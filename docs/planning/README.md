# Planning & Implementation Skills

Six skills that take you from fuzzy requirements to implemented code in an isolated worktree.

## The Pipeline

```
                                    ┌──────────────────────────────────────────────────┐
(optional)          (optional)      │              Core Pipeline                       │
make-scenarios ──→ make-req ──────→ │ make-plan → dry-walkthrough → implement-worktree │
                        ↑           └──────────────────────────────────────────────────┘
                  or use directly
```

**Core pipeline:** `make-plan → dry-walkthrough → implement-worktree`. These three are the minimum for any planned implementation.

**Optional upstream:** `make-scenarios` and `make-req` help when you don't yet know what to build.

**Optional elaboration:** `elaborate-phase` splits large plans into independent sub-plans per phase.

## Skills

### `/make-scenarios` (optional)

Discovers "Actor wants to..." scenarios from a stated perspective. You point at a codebase or component, state your perspective (e.g., "I'm working on the auth module for developers"), and it explores the codebase to surface use cases. Useful for finding requirements you haven't thought of.

### `/make-req` (optional)

Decomposes scenarios (or any task description) into grouped, verifiable requirements (REQ-GRP-NNN). Can also reverse-engineer requirements from an existing codebase. Each requirement is a single testable statement.

### `/make-plan`

Creates implementation plans with architecture diagrams. Uses subagents to explore the codebase, understand existing patterns, and design a solution. Output includes a mermaid diagram (auto-selecting the most appropriate arch lens), tests to write first, ordered implementation steps, and verification criteria.

**Output:** `temp/make-plan/{task_name}_plan_{date}.md`

### `/dry-walkthrough`

Validates a plan by tracing every proposed change against the actual codebase without implementing anything. Catches issues like: files that don't exist, functions with wrong signatures, missing imports, incorrect assumptions about data flow. Fixes the plan directly and reports what changed.

### `/elaborate-phase`

Takes a single phase from a large plan and produces a complete, self-contained implementation plan for just that phase. Includes its own tests, verification, and architecture diagram. Use when a plan is too large to implement in one session.

### `/implement-worktree`

Implements a plan in an isolated git worktree branched from the current branch. Creates the worktree, uses subagents to understand affected systems, implements phase by phase with commits, runs tests, and rebases for squash-and-merge. Never touches the main working directory.

## When to Use What

| Situation | Start with |
|-----------|------------|
| You know exactly what to build | `/make-plan` |
| You have a vague idea or feature area | `/make-scenarios → /make-req → /make-plan` |
| You have a bug to fix | `/investigate → /rectify` (investigation pipeline) |
| Plan is too large for one session | `/elaborate-phase` on individual phases |
| Plan needs external validation | `/review-approach` before implementing |

## Examples

- [Sample make-plan output](examples/make-plan-example.md)
