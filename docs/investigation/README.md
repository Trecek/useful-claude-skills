# Investigation Skills

Three skills for finding, understanding, and solving bugs at the architectural level.

## The Pipeline

```
/investigate → /rectify → /implement-worktree
                  ↑
          (optional: /review-approach)
```

## Skills

### `/investigate`

Deep root cause analysis using parallel subagents. Point it at an error traceback, a bug description, or a codebase question. It spawns subagents to explore affected components, data flow, test gaps, and similar patterns concurrently.

**Output:** `temp/investigate/investigation_{topic}_{date}.md` with sections for root cause, affected components, data flow, test gap analysis, similar patterns, and recommendations.

**When to use:** You have a bug, error, or unexpected behavior and need to understand why it happens before deciding how to fix it.

### `/rectify`

Designs architectural immunity — structural changes that make the bug category impossible, not just the specific instance. Run this after `/investigate`. It identifies test gaps and architectural weaknesses, then produces a plan that addresses root causes.

**Output:** An implementation plan focused on architectural guards (contracts, validation layers, structural constraints) rather than direct bug patches.

**When to use:** After an investigation, when you want a fix that prevents the entire class of bug, not just the one you found.

### `/review-approach`

Researches modern solutions and patterns via web search. Run this on a rectify plan or make-plan to validate that your proposed approach aligns with current best practices. It checks whether better libraries, patterns, or techniques exist for the problem you're solving.

**Output:** Research findings with source URLs, comparison of approaches, and recommendations.

**When to use:** Before implementing a plan, when you want external validation that your approach is sound. Optional but useful for unfamiliar domains.

## Key Design Principle

The pipeline separates **understanding** (investigate) from **solving** (rectify). This prevents the common failure mode of jumping to a fix before understanding the problem. Investigate never modifies code — it only produces analysis. Rectify never investigates — it takes an investigation as input and designs a solution.

## Examples

- [Sample investigation output](examples/investigate-example.md)
