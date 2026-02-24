<a id="top"></a>

# Useful Claude Code Skills

A collection of reusable [Claude Code](https://docs.anthropic.com/en/docs/claude-code) skills for software engineering workflows. I adapted these skills from a specific project to be more generic and project-agnostic.

## What Are Skills?

Claude Code skills are markdown instruction files that live in `.claude/skills/` and teach Claude specific workflows. In this repo, the skills are in that standard location. When invoked (e.g., `/investigate`), Claude follows the skill's methodology — using subagents for parallel exploration, writing structured outputs, and chaining skills together.

[Installation](#installation)

## Skill Catalog

| Category | Skills | What It Does |
|----------|--------|-------------|
| [Architecture Lenses](#architecture-lenses) | 13 | Visualize your codebase from different perspectives using mermaid diagrams |
| [Investigation](#investigation) | 3 | Deep codebase analysis, root cause finding, architectural immunity |
| [Planning & Implementation](#planning--implementation) | 6 | Requirements → plan → validate → implement pipeline |
| [Auditing](#auditing) | 7 | Audit architecture, tests, bug patterns, AI-generated slop, and workflow friction |
| [Documentation](#documentation) | 3 | Keep architecture docs and specs in sync with code |

---

### Architecture Lenses

Visualize your codebase from different architectural perspectives using mermaid diagrams. Each lens answers a specific question about your system.
Use the vscode extension mermaid to see rendered mermaid plots within the markdown files.

| Skill | Lens | Question It Answers |
|-------|------|---------------------|
| [`arch-lens-c4-container`](.claude/skills/arch-lens-c4-container/SKILL.md) | C4 Container | How is it built? |
| [`arch-lens-process-flow`](.claude/skills/arch-lens-process-flow/SKILL.md) | Process Flow | How does it behave? |
| [`arch-lens-data-lineage`](.claude/skills/arch-lens-data-lineage/SKILL.md) | Data Lineage | Where is the data? |
| [`arch-lens-module-dependency`](.claude/skills/arch-lens-module-dependency/SKILL.md) | Module Dependency | How are modules coupled? |
| [`arch-lens-concurrency`](.claude/skills/arch-lens-concurrency/SKILL.md) | Concurrency | How does parallelism work? |
| [`arch-lens-error-resilience`](.claude/skills/arch-lens-error-resilience/SKILL.md) | Error/Resilience | How are failures handled? |
| [`arch-lens-repository-access`](.claude/skills/arch-lens-repository-access/SKILL.md) | Repository Access | How is data accessed? |
| [`arch-lens-operational`](.claude/skills/arch-lens-operational/SKILL.md) | Operational | How is it run and monitored? |
| [`arch-lens-security`](.claude/skills/arch-lens-security/SKILL.md) | Security | Where are the trust boundaries? |
| [`arch-lens-development`](.claude/skills/arch-lens-development/SKILL.md) | Development | How is it built and tested? |
| [`arch-lens-scenarios`](.claude/skills/arch-lens-scenarios/SKILL.md) | Scenarios | Do the components work together? |
| [`arch-lens-state-lifecycle`](.claude/skills/arch-lens-state-lifecycle/SKILL.md) | State Lifecycle | How is state corruption prevented? |
| [`arch-lens-deployment`](.claude/skills/arch-lens-deployment/SKILL.md) | Deployment | Where does it run? |

Use [`make-arch-diag`](.claude/skills/make-arch-diag/SKILL.md) to select the right lens interactively. Use [`verify-diag`](.claude/skills/verify-diag/SKILL.md) to validate a diagram's accuracy against the codebase. 
The make-plan and rectify skills will automatically choose the most appropriate lens for visualizing the proposed changes to codebase.

See [examples generated against UMI-tools](docs/arch-lens/examples/umi-tools/) ([process-flow](docs/arch-lens/examples/umi-tools/process-flow.md), [data-lineage](docs/arch-lens/examples/umi-tools/data-lineage.md), [module-dependency](docs/arch-lens/examples/umi-tools/module-dependency.md))

---

### Investigation

Skills for investigating bugs, making plans to address bugs, and researching best practices & patterns. The core flow is investigate a problem, then design a solution that solves the architecture rather than the bug. `review-approach` is an optional step you can run on a rectify plan or make-plan to research what modern solutions exist.

| Skill | Purpose | Trigger |
|-------|---------|---------|
| [`investigate`](.claude/skills/investigate/SKILL.md) | Root cause analysis with parallel subagents | `/investigate`, paste an error traceback |
| [`rectify`](.claude/skills/rectify/SKILL.md) | Devise architectural immunity plans (not bandaid fixes) | `/rectify` after an investigation |
| [`review-approach`](.claude/skills/review-approach/SKILL.md) | Research modern solutions via web search (optional) | `/review-approach` on any plan |

---

### Planning & Implementation

The core pipeline is **make-plan → dry-walkthrough → implement-worktree**. The make-scenarios & make-req skills are optional exploratory steps for when you're not yet sure what to build.
elaborate-phase can be used when your plan is too large to implement in one go. It will make an independent plan for each phase.

```
                                    ┌──────────────────────────────────────────────────┐
(optional)          (optional)      │              Core Pipeline                       │
make-scenarios ──→ make-req ──────→ │ make-plan → dry-walkthrough → implement-worktree │
                        ↑           └──────────────────────────────────────────────────┘
                  or use directly
```

**`make-scenarios`** helps when you don't yet know how to frame a problem or feature. You state a perspective (e.g., "I'm working on the authentication module for developers") and the skill explores the codebase to surface scenarios like "Developer doesn't want to get logged out mid-task." Those scenarios then guide requirement writing.

**`make-req`** takes scenarios (or any task description) and decomposes them into grouped, verifiable requirements. It can also be pointed at an entire repo to reverse-engineer the requirements that would be needed to generate the project from scratch.

Both are useful for refining what you need before planning. Neither is required — you can go straight to `/make-plan` if you already know what to build.

| Skill | Purpose | Trigger |
|-------|---------|---------|
| [`make-scenarios`](.claude/skills/make-scenarios/SKILL.md) | Discover "Actor wants to..." scenarios from a stated perspective | `/make-scenarios` |
| [`make-req`](.claude/skills/make-req/SKILL.md) | Decompose into grouped, verifiable requirements (REQ-GRP-NNN) | `/make-req` |
| [`make-plan`](.claude/skills/make-plan/SKILL.md) | Create implementation plans with arch lens diagrams | `/make-plan` |
| [`elaborate-phase`](.claude/skills/elaborate-phase/SKILL.md) | Elaborate a single phase into a self-contained implementation plan | `/elaborate-phase` |
| [`dry-walkthrough`](.claude/skills/dry-walkthrough/SKILL.md) | Validate a plan by tracing every change against the codebase | `/dry-walkthrough` |
| [`implement-worktree`](.claude/skills/implement-worktree/SKILL.md) | Implement a plan in an isolated git worktree | `/implement-worktree` |

---

### Auditing

Audit codebases for architectural issues, test quality, bug patterns, and AI-generated slop.
These should be tailored to your specific project.

| Skill | Purpose | Trigger |
|-------|---------|---------|
| [`audit-arch`](.claude/skills/audit-arch/SKILL.md) | Audit architecture against standards and practices | `/audit-arch` |
| [`audit-bugs`](.claude/skills/audit-bugs/SKILL.md) | Mine historical bug patterns from claude code project conversation logs | `/audit-bugs` |
| [`audit-defense-standards`](.claude/skills/audit-defense-standards/SKILL.md) | Audit codebase against defense standards from bug patterns | `/audit-defense-standards` |
| [`audit-tests`](.claude/skills/audit-tests/SKILL.md) | Find useless tests, over-mocking, weak assertions | `/audit-tests` |
| [`design-guards`](.claude/skills/design-guards/SKILL.md) | Design architectural guards for identified bug patterns | `/design-guards` |
| [`id-slop`](.claude/skills/id-slop/SKILL.md) | Find AI-generated code slop (phase comments, dead code, compat hacks) | `/id-slop` |
| [`audit-friction`](.claude/skills/audit-friction/SKILL.md) | Mine session logs for repeated failures, stuck loops, and workflow friction | `/audit-friction` |

---

### Documentation

Keep architecture docs and specifications in sync with implementation.
These should be tailored to your specific project.

| Skill | Purpose | Trigger |
|-------|---------|---------|
| [`update-architecture`](.claude/skills/update-architecture/SKILL.md) | Rebuild architecture docs from code understanding (8-file structure) | `/update-architecture` |
| [`update-specs`](.claude/skills/update-specs/SKILL.md) | Maintain append-only functional specifications (SPEC-NNN) | `/update-specs` |
| [`mermaid`](.claude/skills/mermaid/SKILL.md) | Standard mermaid diagram styling and conventions | `/mermaid` |

---

## How Skills Work Together

Many skills are designed to chain. Here are common workflows:

**Bug discovered → architectural fix:**
```
/investigate → /rectify → /dry-walkthrough → /implement-worktree
                  ↑
          (optional: /review-approach on the rectify plan)
```

**New feature → shipped:**
```
/make-plan → /dry-walkthrough → /implement-worktree
     ↑
(optional: /make-scenarios → /make-req to clarify what to build)
(optional: /review-approach on the plan)
(optional: /elaborate-phase on the plan to split into smaller independent plans)
```

**Codebase health check → fixes:**
```
/audit-arch + /audit-tests + /id-slop → /make-plan → /dry-walkthrough → /implement-worktree
```

## Installation

To install for a specific project, open a terminal in your project folder and run:

```bash
git clone https://github.com/Trecek/useful-claude-skills.git /tmp/claude-skills && \
  mkdir -p .claude/skills && \
  cp -r /tmp/claude-skills/.claude/skills/* .claude/skills/ && \
  rm -rf /tmp/claude-skills
```

To install globally (available in all projects):

```bash
git clone https://github.com/Trecek/useful-claude-skills.git /tmp/claude-skills && \
  mkdir -p ~/.claude/skills && \
  cp -r /tmp/claude-skills/.claude/skills/* ~/.claude/skills/ && \
  rm -rf /tmp/claude-skills
```

To install a single skill (e.g., `investigate`):

```bash
git clone --depth 1 https://github.com/Trecek/useful-claude-skills.git /tmp/claude-skills && \
  mkdir -p .claude/skills/investigate && \
  cp /tmp/claude-skills/.claude/skills/investigate/SKILL.md .claude/skills/investigate/ && \
  rm -rf /tmp/claude-skills
```

[Back to top](#top)
