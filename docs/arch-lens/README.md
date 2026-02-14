# Architecture Lens Skills

Mermaid diagrams that visualize your codebase from 13 different architectural perspectives. Each lens answers a specific question about your system.

## How Lenses Work

A lens is a focused perspective on your architecture. You don't use all 13 on every project — you pick the lens that answers the question you're asking right now.

**Scoping:** Lenses are scoped to the subsystem being changed, not the whole project. A diagram for "add retry logic to the API client" should show the API client's error handling path, not the entire application.

## Lens Selection Guide

| Your plan involves... | Use this lens | Question it answers |
|----------------------|---------------|---------------------|
| Adding/restructuring components | C4 Container | How is it built? |
| Changing runtime behavior | Process Flow | How does it behave? |
| Modifying data pipelines | Data Lineage | Where is the data? |
| Refactoring module structure | Module Dependency | How are modules coupled? |
| Adding parallelism/threading | Concurrency | How does parallelism work? |
| Changing error handling | Error/Resilience | How are failures handled? |
| Modifying data access patterns | Repository Access | How is data accessed? |
| Changing CLI/ops workflows | Operational | How is it run and monitored? |
| Adding auth/validation | Security | Where are the trust boundaries? |
| Changing build/test setup | Development | How is it built and tested? |
| Validating component cooperation | Scenarios | Do the components work together? |
| Modifying state management | State Lifecycle | How is state corruption prevented? |
| Changing deployment/infra | Deployment | Where does it run? |

## How Lenses Chain

```
make-plan / rectify
    → auto-selects appropriate arch lens
    → lens loads mermaid skill for styling
    → verify-diag validates against codebase
```

The `make-plan` and `rectify` skills automatically choose the most appropriate lens based on what the plan modifies. You can also use `make-arch-diag` to interactively select a lens.

## Shared Styling

All lenses use the same mermaid `classDef` color palette defined in the [mermaid skill](../../.claude/skills/mermaid/SKILL.md). This ensures visual consistency across diagrams.

## Examples

### UMI-tools (Python bioinformatics CLI)

| Lens | Example |
|------|---------|
| C4 Container | [c4-container.md](examples/umi-tools/c4-container.md) |
| Process Flow | [process-flow.md](examples/umi-tools/process-flow.md) |
| Data Lineage | [data-lineage.md](examples/umi-tools/data-lineage.md) |
| Module Dependency | [module-dependency.md](examples/umi-tools/module-dependency.md) |
| Concurrency | [concurrency.md](examples/umi-tools/concurrency.md) |
| Error/Resilience | [error-resilience.md](examples/umi-tools/error-resilience.md) |
| Operational | [operational.md](examples/umi-tools/operational.md) |
| Security | [security.md](examples/umi-tools/security.md) |

### STAR Aligner (C++ genomic aligner)

| Lens | Example |
|------|---------|
| C4 Container | [c4-container.md](examples/star/c4-container.md) |
| Concurrency | [concurrency.md](examples/star/concurrency.md) |
| Development | [development.md](examples/star/development.md) |

### Nextflow (Groovy/Java workflow engine)

| Lens | Example |
|------|---------|
| Deployment | [deployment.md](examples/nextflow/deployment.md) |
| Process Flow | [process-flow.md](examples/nextflow/process-flow.md) |
| Data Lineage | [data-lineage.md](examples/nextflow/data-lineage.md) |
| Scenarios | [scenarios.md](examples/nextflow/scenarios.md) |

### Scanpy (Python single-cell analysis)

| Lens | Example |
|------|---------|
| Data Lineage | [data-lineage.md](examples/scanpy/data-lineage.md) |
| Module Dependency | [module-dependency.md](examples/scanpy/module-dependency.md) |
| State Lifecycle | [state-lifecycle.md](examples/scanpy/state-lifecycle.md) |
| Repository Access | [repository-access.md](examples/scanpy/repository-access.md) |
