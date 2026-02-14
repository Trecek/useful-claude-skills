# Useful Claude Code Skills

A collection of reusable [Claude Code](https://docs.anthropic.com/en/docs/claude-code) skills for software engineering workflows. Each skill is project-agnostic and can be dropped into any codebase.

## What Are Skills?

Claude Code skills are markdown instruction files that live in `.claude/skills/` and teach Claude specific workflows. When invoked (e.g., `/investigate`), Claude follows the skill's methodology — using subagents for parallel exploration, writing structured outputs, and chaining skills together.

## Installation

Skills live in a `.claude/skills/` folder. You can install them **per-project** (only available in that project) or **globally** (available in all your projects).

### Quick Install (easiest)

Open a terminal in your project folder and run:

```bash
# Install ALL skills into your current project
git clone https://github.com/Trecek/useful-claude-skills.git /tmp/claude-skills && \
  mkdir -p .claude/skills && \
  cp -r /tmp/claude-skills/skills/* .claude/skills/ && \
  rm -rf /tmp/claude-skills
```

That's it. The skills are now available when you use Claude Code in this project.

### Install Globally (all projects)

To make skills available everywhere, install to `~/.claude/skills/`:

```bash
git clone https://github.com/Trecek/useful-claude-skills.git /tmp/claude-skills && \
  mkdir -p ~/.claude/skills && \
  cp -r /tmp/claude-skills/skills/* ~/.claude/skills/ && \
  rm -rf /tmp/claude-skills
```

### Install a Single Skill

If you only want one skill (e.g., `investigate`):

```bash
git clone --depth 1 https://github.com/Trecek/useful-claude-skills.git /tmp/claude-skills && \
  mkdir -p .claude/skills/investigate && \
  cp /tmp/claude-skills/skills/investigate/SKILL.md .claude/skills/investigate/ && \
  rm -rf /tmp/claude-skills
```

### Download Without Git

If you don't have git installed, download the zip from GitHub:

1. Go to https://github.com/Trecek/useful-claude-skills
2. Click the green **Code** button, then **Download ZIP**
3. Unzip the file
4. Copy the contents of the `skills/` folder into your project's `.claude/skills/` folder (or `~/.claude/skills/` for global)

### What Gets Installed

Each skill is a single `SKILL.md` file inside its own folder. The structure looks like:

```
.claude/skills/
├── investigate/
│   └── SKILL.md
├── make-plan/
│   └── SKILL.md
├── audit-tests/
│   └── SKILL.md
└── ...
```

Skills are self-contained — each `SKILL.md` has YAML frontmatter for hooks and a markdown body with instructions. No dependencies, no build steps.

## Skill Catalog

### Architecture Lenses (13 skills)

Visualize your codebase from different architectural perspectives using mermaid diagrams. Each lens answers a specific question about your system.

| Skill | Lens | Question It Answers |
|-------|------|---------------------|
| [`arch-lens-c4-container`](skills/arch-lens-c4-container/SKILL.md) | C4 Container | How is it built? |
| [`arch-lens-process-flow`](skills/arch-lens-process-flow/SKILL.md) | Process Flow | How does it behave? |
| [`arch-lens-data-lineage`](skills/arch-lens-data-lineage/SKILL.md) | Data Lineage | Where is the data? |
| [`arch-lens-module-dependency`](skills/arch-lens-module-dependency/SKILL.md) | Module Dependency | How are modules coupled? |
| [`arch-lens-concurrency`](skills/arch-lens-concurrency/SKILL.md) | Concurrency | How does parallelism work? |
| [`arch-lens-error-resilience`](skills/arch-lens-error-resilience/SKILL.md) | Error/Resilience | How are failures handled? |
| [`arch-lens-repository-access`](skills/arch-lens-repository-access/SKILL.md) | Repository Access | How is data accessed? |
| [`arch-lens-operational`](skills/arch-lens-operational/SKILL.md) | Operational | How is it run and monitored? |
| [`arch-lens-security`](skills/arch-lens-security/SKILL.md) | Security | Where are the trust boundaries? |
| [`arch-lens-development`](skills/arch-lens-development/SKILL.md) | Development | How is it built and tested? |
| [`arch-lens-scenarios`](skills/arch-lens-scenarios/SKILL.md) | Scenarios | Do the components work together? |
| [`arch-lens-state-lifecycle`](skills/arch-lens-state-lifecycle/SKILL.md) | State Lifecycle | How is state corruption prevented? |
| [`arch-lens-deployment`](skills/arch-lens-deployment/SKILL.md) | Deployment | Where does it run? |

Use [`make-arch-diag`](skills/make-arch-diag/SKILL.md) to select the right lens interactively.

[More details and examples](docs/arch-lens/)

---

### Investigation (3 skills)

Deep codebase analysis without making changes. The core flow is investigate a problem, then devise architectural immunity. `review-approach` is an optional step you can run on a rectify plan or make-plan to research what modern solutions exist before committing to a direction.

| Skill | Purpose | Trigger |
|-------|---------|---------|
| [`investigate`](skills/investigate/SKILL.md) | Root cause analysis with parallel subagents | `/investigate`, paste an error traceback |
| [`rectify`](skills/rectify/SKILL.md) | Devise architectural immunity plans (not bandaid fixes) | `/rectify` after an investigation |
| [`review-approach`](skills/review-approach/SKILL.md) | Research modern solutions via web search (optional) | `/review-approach` on any plan |

[More details](docs/investigation/)

---

### Planning & Implementation (6 skills)

The core pipeline is **make-plan → elaborate-phase → dry-walkthrough → implement-worktree**. The first two skills are optional exploratory steps for when you're not yet sure what to build.

```
                                    ┌─────────────────────────────────────────────────────────────────────┐
(optional)          (optional)      │              Core Pipeline                                          │
make-scenarios ──→ make-req ──────→ │ make-plan → elaborate-phase → dry-walkthrough → implement-worktree  │
                        ↑           └─────────────────────────────────────────────────────────────────────┘
                  or use directly
```

**`make-scenarios`** is a Jeopardy-style discovery step — you know what you *want* but not what to *implement*. You state a perspective (e.g., "as a developer, I want a modular codebase with well-defined abstraction layers so I never need to touch unrelated components") and the skill explores the codebase to surface scenarios that would address that desire. It can also identify all the scenarios a developer, user, or operator might have for a specific component.

**`make-req`** takes scenarios (or any task description) and decomposes them into grouped, verifiable requirements. It turns desires into concrete acceptance criteria.

Both are useful for refining what you need before planning. Neither is required — you can go straight to `/make-plan` if you already know what to build.

| Skill | Purpose | Trigger |
|-------|---------|---------|
| [`make-scenarios`](skills/make-scenarios/SKILL.md) | Discover "Actor wants to..." scenarios from a stated perspective | `/make-scenarios` |
| [`make-req`](skills/make-req/SKILL.md) | Decompose into grouped, verifiable requirements (REQ-GRP-NNN) | `/make-req` |
| [`make-plan`](skills/make-plan/SKILL.md) | Create implementation plans with arch lens diagrams | `/make-plan` |
| [`elaborate-phase`](skills/elaborate-phase/SKILL.md) | Elaborate a single phase into a self-contained implementation plan | `/elaborate-phase` |
| [`dry-walkthrough`](skills/dry-walkthrough/SKILL.md) | Validate a plan by tracing every change against the codebase | `/dry-walkthrough` |
| [`implement-worktree`](skills/implement-worktree/SKILL.md) | Implement a plan in an isolated git worktree | `/implement-worktree` |

[More details](docs/planning/)

---

### Auditing (6 skills)

Audit codebases for architectural issues, test quality, bug patterns, and AI-generated slop.

| Skill | Purpose | Trigger |
|-------|---------|---------|
| [`audit-arch`](skills/audit-arch/SKILL.md) | Audit architecture against standards and practices | `/audit-arch` |
| [`audit-bugs`](skills/audit-bugs/SKILL.md) | Mine historical bug patterns from investigation logs | `/audit-bugs` |
| [`audit-defense-standards`](skills/audit-defense-standards/SKILL.md) | Audit codebase against defense standards from bug patterns | `/audit-defense-standards` |
| [`audit-tests`](skills/audit-tests/SKILL.md) | Find useless tests, over-mocking, weak assertions | `/audit-tests` |
| [`design-guards`](skills/design-guards/SKILL.md) | Design architectural guards for identified bug patterns | `/design-guards` |
| [`id-slop`](skills/id-slop/SKILL.md) | Find AI-generated code slop (phase comments, dead code, compat hacks) | `/id-slop` |

[More details](docs/auditing/)

---

### Documentation (3 skills)

Keep architecture docs and specifications in sync with implementation.

| Skill | Purpose | Trigger |
|-------|---------|---------|
| [`update-architecture`](skills/update-architecture/SKILL.md) | Rebuild architecture docs from code understanding (8-file structure) | `/update-architecture` |
| [`update-specs`](skills/update-specs/SKILL.md) | Maintain append-only functional specifications (SPEC-NNN) | `/update-specs` |
| [`mermaid`](skills/mermaid/SKILL.md) | Standard mermaid diagram styling and conventions | `/mermaid` |

[More details](docs/documentation/)

## How Skills Work Together

Many skills are designed to chain. Here are common workflows:

**Bug discovered → architectural fix:**
```
/investigate → /rectify → /make-plan → /dry-walkthrough → /implement-worktree
                  ↑
          (optional: /review-approach on the rectify plan)
```

**New feature → shipped:**
```
/make-plan → /elaborate-phase → /dry-walkthrough → /implement-worktree
     ↑
(optional: /make-scenarios → /make-req to clarify what to build)
(optional: /review-approach on the plan)
```

**Codebase health check → fixes:**
```
/audit-arch + /audit-tests + /id-slop → /make-plan → /dry-walkthrough → /implement-worktree
```

**Post-implementation:**
```
/update-architecture → /update-specs
```

## Customization

These skills are generic starting points. To tailor them to your project:

- **Test commands**: Replace "run the project's test suite" with your actual command (e.g., `pytest`, `npm test`, `cargo test`)
- **Output directories**: Skills write to `temp/{skill-name}/` by default — adjust if your project uses a different temp location
- **Lint/format checks**: Add your specific linter commands to verification steps
- **Architecture docs path**: Update `update-architecture` with your docs directory structure

## License

MIT
