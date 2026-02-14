# Documentation Skills

Three skills for keeping documentation in sync with code.

## Skills

### `/update-architecture`

Rebuilds architecture documentation from direct code understanding. Maintains an 8-file structure:

| File | Purpose |
|------|---------|
| `README.md` | Navigation hub and quick reference |
| `overview.md` | System diagrams and design principles |
| `components.md` | Component details and dependencies |
| `data-models.md` | State schemas and database models |
| `workflows.md` | Execution flows and state machines |
| `patterns.md` | Design patterns, ID formats, conventions |
| `configuration.md` | CLI flags, env vars, settings |
| `directory-structure.md` | Code navigation and module layout |

Each file is updated selectively — only files affected by recent changes get rewritten. The skill reads the codebase directly rather than relying on stale docs.

### `/update-specs`

Maintains append-only functional specifications using `SPEC-NNN` format. Each spec is a single testable requirement statement. Specs are never deleted — obsolete ones are marked with `Status: Warning`. Numbers never change, ensuring stable references.

### `/mermaid`

Shared mermaid diagram styling and conventions used by all diagram-producing skills. Defines the standard `classDef` color palette, node label format, init configuration, and layout rules. Loaded automatically by arch-lens skills and make-plan.

## Workflow

Run `/update-architecture` and `/update-specs` after implementing changes to keep docs current:

```
/implement-worktree → /update-architecture → /update-specs
```

## Examples

- [Sample architecture doc structure](examples/update-architecture-example.md)
