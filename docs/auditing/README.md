# Auditing Skills

Six skills for assessing codebase health. Split into two modes: proactive (find issues before they cause bugs) and reactive (learn from bugs that already happened).

## Proactive Auditing

```
/audit-arch + /audit-tests + /id-slop → findings feed /make-plan
```

### `/audit-arch`

Audits codebase against architectural standards and practices. Spawns parallel subagents to examine multiple architectural aspects (layering violations, dependency direction, naming conventions, pattern consistency). Generates a structured report with findings and severity.

### `/audit-tests`

Finds test quality issues across 9 categories: useless tests, redundant/consolidatable tests, over-mocked tests, weak assertions, misleading tests, stale tests, fixture issues, misclassified tests, and oversized test files. Each finding includes file path, line range, severity, explanation, and recommended action.

**Output:** `temp/audit-tests/test_audit_{date}.md`

### `/id-slop`

Identifies AI-generated code slop: useless comments (phase markers, section dividers), backward compatibility hacks, deprecation notices for code that was never public, dead code paths, and unnecessary abstractions. Generates a removal plan validated by dry walkthrough.

## Reactive Auditing

```
/audit-bugs → /design-guards → /audit-defense-standards
```

### `/audit-bugs`

Mines historical bug patterns from Claude Code project conversation logs. Identifies recurring root causes, categorizes failure modes, and surfaces architectural gaps that produce repeated bugs. Use when you want to understand what keeps breaking and why.

### `/design-guards`

Takes an audit-bugs report and designs architectural guards (tests, contracts, structural changes) that provide immunity to each identified pattern. Guards are structural — they make the bug category impossible rather than catching individual instances.

### `/audit-defense-standards`

Audits the codebase against defense standards derived from historical bug patterns. Standards accumulate over time as new patterns are discovered via `/audit-bugs` and `/design-guards`. This skill needs project-specific customization to define what standards apply.

## Typical Workflow

1. Run `/audit-arch`, `/audit-tests`, and `/id-slop` periodically for proactive health checks
2. When bugs cluster, run `/audit-bugs` to find patterns
3. Run `/design-guards` on the bug report to create structural defenses
4. Use findings from any audit as input to `/make-plan` for implementation

## Examples

- [Sample audit-tests output](examples/audit-tests-example.md)
