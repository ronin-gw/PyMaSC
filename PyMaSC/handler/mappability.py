import logging
import os
import json

from PyMaSC.core.mappability import BWFeederWithMappableRegionSum

logger = logging.getLogger(__name__)


class BWIOError(IOError):
    pass


class JSONIOError(IOError):
    pass


class NeedUpdate(Exception):
    pass


class MappabilityHandler(BWFeederWithMappableRegionSum):
    def __init__(self, path, max_shift=0, map_path=None):
        if not os.path.isfile(path):
            logger.critical("Failed to open '{}': no such file.".format(path))
            raise BWIOError
        elif not os.access(path, os.R_OK):
            logger.critical("Failed to open '{}': file is unreadable.".format(path))
            raise BWIOError

        super(MappabilityHandler, self).__init__(path, max_shift)
        self.need_save_stats = True

        if map_path:
            self.map_path = map_path
        else:
            self.map_path = os.path.splitext(path)[0] + "_mappability.json"

        if not os.path.exists(self.map_path):
            dirname = os.path.dirname(self.map_path)
            dirname = dirname if dirname else '.'
            if not os.access(dirname, os.W_OK):
                logger.critical("Directory is not writable: '{}'".format(dirname))
                raise JSONIOError
        elif not os.path.isfile(self.map_path):
            logger.critical("Specified path is not file: '{}'".format(self.map_path))
            raise JSONIOError
        else:
            if not os.access(self.map_path, os.R_OK):
                logger.error("Failed to read '{}'".format(self.map_path))
            else:
                logger.info("Read mappability stats from '{}'".format(self.map_path))
                try:
                    self._load_mappability_stats()
                except IOError as e:
                    logger.error("Failed to read '{}'".format(self.map_path))
                    logger.error("[Errno {}] {}".format(e.errno, e.message))
                except (TypeError, OverflowError, ValueError, KeyError, IndexError) as e:
                    logger.error("Failed to load json file: '{}'".format(self.map_path))
                except NeedUpdate:
                    pass

                if self.need_save_stats:
                    if not os.access(self.map_path, os.W_OK):
                        logger.critical("Failed to overwrite '{}'".format(self.map_path))
                        raise JSONIOError
                    else:
                        logger.warning("Existing file '{}' will be overwritten.".format(self.map_path))
                    logger.info("Calcurate mappable length with max shift size {}.".format(max_shift))

    def _load_mappability_stats(self):
        with open(self.map_path) as f:
            stats = json.load(f)

        for k in ("max_shift", "__whole__", "references"):
            if k not in stats:
                logger.error("Mandatory key '{}' not found.".format(k))
                raise KeyError(k)

        if stats["max_shift"] < self.max_shift:
            logger.info("Specified shift length longer than former analysis. The stats will be updated.")
            raise NeedUpdate

        if stats["max_shift"] != len(stats["__whole__"]) - 1:
            logger.error("Max shift length for whole genome unmatched.")
            raise IndexError

        for ref in self.chromsizes:
            if ref not in stats["references"]:
                logger.error("Reference '{}' not found.".format(ref))
                raise KeyError(ref)

            if stats["max_shift"] != len(stats["references"][ref]) - 1:
                logger.error("Max shift length for 'ref' unmatched.".format(ref))
                raise IndexError

        self.mappable_len = stats["__whole__"][:self.max_shift + 1]
        self.chrom2mappable_len = {ref: b[:self.max_shift + 1] for ref, b in stats["references"].items()}
        self.chrom2is_called = {ref: True for ref in self.chromsizes}
        self.is_called = True
        self.need_save_stats = False

    def save_mappability_stats(self):
        if not self.need_save_stats:
            return logger.info("Update mappability stats is not required.")

        self.calc_mappability()

        logger.info("Save mappable length to '{}'".format(self.map_path))
        try:
            with open(self.map_path, 'w') as f:
                json.dump({
                    "max_shift": self.max_shift,
                    "__whole__": self.mappable_len,
                    "references": self.chrom2mappable_len
                }, f, indent=4, sort_keys=True)
        except IOError as e:
            logger.error("Faild to output: {}\n[Errno {}] {}".format(
                         e.filename, e.errno, e.message))

        self.need_save_stats = False
