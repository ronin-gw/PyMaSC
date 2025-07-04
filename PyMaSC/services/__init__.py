"""Service layer implementations for PyMaSC.

This package provides service layer abstractions that separate
different concerns into focused, testable components:

- calculation: Pure calculation logic (domain layer)
- io: Input/output operations (infrastructure layer)
- workflow: Process orchestration (application layer)

The service layer architecture improves:
- Testability through dependency injection
- Maintainability through single responsibility
- Flexibility through interface-based design
"""