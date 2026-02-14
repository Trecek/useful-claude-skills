# **Project Development Guidelines**

This document provides mandatory instructions for AI-assisted development in this repository. All contributions must be consistent, high-quality, and aligned with these standards.

## **1. Core Project Goal**

<!-- Replace with your project's objective -->
_Describe the primary objective of this project here._

## **2. General Principles**

  * **Follow the Task Description**: Your primary source of truth is the issue, ticket, or work package description provided for your assignment. It contains the complete scope and requirements for your work.
  * **Adhere to Task Scope**: Strictly adhere to the scope of the assigned task. Do not work on unassigned features, unrelated refactoring, or bug fixes outside the current assignment.
  * **Implement Faithfully**: Produce a functionally correct implementation based on the task requirements. Do not add new features or architectural changes unless explicitly requested.
  * **Adhere to Project Standards**: Write clean, maintainable code following the established conventions and architectural patterns of this project.

## **3. Critical Rules - DO NOT VIOLATE**

These rules are essential for maintaining project structure, preventing bugs, and ensuring a clean codebase.

### **3.1. Code and Implementation**

  * **Do Not Oversimplify**: Implement logic with its required complexity. Do not take shortcuts that compromise correctness or violate requirements.
  * **Respect the Existing Architecture**: Build upon the established project structure, modules, and design patterns. Do not introduce new architectural patterns or file structures without explicit instruction. Always understand existing architecture before implementing a new feature.
  * **Address the Root Cause**: When your code fails a test or causes a bug, debug the implementation to find and fix the root cause. Avoid hardcoded workarounds that only solve for specific inputs.
  * **No Backward Compatibility Hacks**: Do NOT leave comments about dead or deprecated code. No backward compatibility or deprecated code should be kept. Dead code should always be removed.
  * **Avoid Redundancy**: Do not duplicate logic, utilities, or test code.
  * **Use Current Package Versions**: When adding new dependencies, use web search to find and use the current stable version.

### **3.2. File System**

  * **Temporary Files:** All temporary files for debugging, analysis, or testing hypotheses **must** be created in the project's `temp/` directory.
  * **Do Not Add Root Files**: Never create new files or scripts in the repository root unless the task explicitly requires it.
  * **Never commit unless told to do so**

## **4. Testing Guidelines**

### **4.1. Running Tests**

<!-- Replace with your project's test commands -->
```bash
# Example: pytest, npm test, cargo test, etc.
```

### **4.2. Test Requirements**

  * **Always run tests at end of task**
  * **Fix failing tests immediately**
  * **Add tests for new features**: Every new feature should have corresponding tests.
  * **Follow existing test patterns**: New tests should follow existing testing framework and architecture. Test redundancy must be avoided.

## **5. Pre-commit Hooks**

Install hooks after cloning: `pre-commit install`

Hooks run automatically on commit. To run manually: `pre-commit run --all-files`

<!-- List your project's specific hooks here -->

## **6. Architecture**

<!-- Replace with your project's architecture overview and directory structure -->
