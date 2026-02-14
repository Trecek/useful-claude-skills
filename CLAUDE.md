# **Project Development Guidelines**

This document provides mandatory instructions for AI-assisted development in this repository. All contributions must be consistent, high-quality, and aligned with these standards.

## **1. Core Project Goal**

A collection of reusable Claude Code skills for software engineering workflows. Skills are markdown instruction files in `.claude/skills/` that teach Claude specific methodologies — investigation, planning, auditing, documentation, and architecture visualization. They chain together (e.g., `investigate → rectify → implement-worktree`) and use subagents for parallel exploration.

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

This is a skills-only repository with no executable code and no test suite. Verification is done by reviewing skill outputs against the format defined in each `SKILL.md`.

## **5. Pre-commit Hooks**

Install hooks after cloning: `pre-commit install`

Hooks run automatically on commit. To run manually: `pre-commit run --all-files`

<!-- List your project's specific hooks here -->

## **6. Architecture**

```
.claude/skills/           # All skill definitions (SKILL.md files)
  arch-lens-*/            # 13 architecture lens skills
  investigate/            # Root cause analysis
  rectify/                # Architectural immunity plans
  review-approach/        # External solution research
  make-plan/              # Implementation planning
  make-scenarios/         # Scenario discovery
  make-req/               # Requirement decomposition
  elaborate-phase/        # Phase elaboration
  dry-walkthrough/        # Plan validation
  implement-worktree/     # Isolated implementation
  audit-arch/             # Architecture audit
  audit-bugs/             # Bug pattern mining
  audit-tests/            # Test quality audit
  audit-defense-standards/ # Defense standard audit
  design-guards/          # Guard design from bug patterns
  id-slop/                # AI slop identification
  update-architecture/    # Architecture doc maintenance
  update-specs/           # Specification maintenance
  mermaid/                # Shared mermaid styling
  make-arch-diag/         # Interactive lens selection
  verify-diag/            # Diagram verification
docs/                     # Category documentation and examples
  arch-lens/              # Lens overview + examples by tool
  investigation/          # Investigation pipeline docs
  planning/               # Planning pipeline docs
  auditing/               # Auditing pipeline docs
  documentation/          # Documentation skills docs
temp/                     # Temporary/working files (gitignored)
```
