# Architecture Documentation Structure

> **Example Note:** This is a sample output for demonstration purposes.
> Actual outputs will reflect your specific codebase and issues.

## Overview

The `/update-architecture` skill maintains a comprehensive 8-file architecture documentation system. Each file serves a specific purpose and should be updated at different points in the development lifecycle.

## Documentation Files

| File | Purpose | Update When |
|------|---------|-------------|
| `OVERVIEW.md` | High-level system architecture and design philosophy | Major architectural changes or new subsystems added |
| `COMPONENTS.md` | Detailed component catalog with responsibilities | New components added or significant refactoring |
| `DATA_FLOW.md` | How data moves through the system | New data pipelines or integration points |
| `DEPLOYMENT.md` | Infrastructure, deployment, and operational concerns | Infrastructure changes, new environments, or deployment process updates |
| `DECISIONS.md` | Architectural Decision Records (ADRs) | Significant technical decisions requiring documentation |
| `DEPENDENCIES.md` | External dependencies and integration points | New third-party services, libraries, or API integrations |
| `SECURITY.md` | Security architecture and threat model | Security features, authentication changes, or compliance requirements |
| `TESTING.md` | Testing strategy and quality practices | Test architecture changes or new testing approaches |

---

## Example File Contents

### OVERVIEW.md

```markdown
# System Architecture Overview

## System Purpose
Multi-tenant SaaS platform for batch data processing with real-time monitoring.

## Core Design Principles
1. **Separation of Concerns:** API layer, business logic, and data persistence are isolated
2. **Fail-Safe Processing:** All batch operations are transactional with automatic rollback
3. **Observable by Default:** Every component emits structured logs and metrics

## High-Level Architecture
[Mermaid diagram showing: API Gateway -> Application Layer -> Service Layer -> Data Layer]
```

---

### COMPONENTS.md

```markdown
# System Components

## API Layer
- **Component:** `src/api/` - REST API endpoints using FastAPI
- **Responsibility:** Request validation, authentication, rate limiting
- **Key Dependencies:** FastAPI, Pydantic, JWT authentication middleware

## Batch Processing Engine
- **Component:** `src/batch/processor.py` - Core batch processing orchestration
- **Responsibility:** Manages batch lifecycle from upload through completion
- **Key Dependencies:** Celery for async tasks, Redis for task queue
```

---

### DATA_FLOW.md

```markdown
# Data Flow Documentation

## Batch Processing Flow
1. Client uploads CSV file via POST /api/v1/batches
2. File validator checks encoding, size, and schema
3. Job queued to Celery with batch_id
4. Worker processes records in transactions of 1000 rows
5. Results written to PostgreSQL, metrics to Prometheus
6. Completion notification sent via webhook

## Authentication Flow
1. User submits credentials to /auth/login
2. Credentials validated against PostgreSQL users table
3. JWT token generated with 1-hour expiration
4. Token returned to client, stored in Redis with user session
```

---

### DEPLOYMENT.md

```markdown
# Deployment Architecture

## Infrastructure
- **Platform:** AWS (primary), Docker containers via ECS
- **Database:** PostgreSQL 14 on RDS with Multi-AZ deployment
- **Cache:** Redis 7 on ElastiCache for sessions and task queue
- **Object Storage:** S3 for uploaded batch files and reports

## Environments
- **Production:** `prod.example.com` - Auto-scaling 4-20 instances
- **Staging:** `staging.example.com` - Fixed 2 instances, mirrors production config
- **Development:** Local Docker Compose stack
```

---

### DECISIONS.md

```markdown
# Architectural Decision Records

## ADR-001: Use PostgreSQL for Primary Database
**Date:** 2025-11-15
**Status:** Accepted
**Context:** Need ACID transactions for batch processing integrity
**Decision:** PostgreSQL chosen over MongoDB for strong consistency guarantees
**Consequences:** Excellent data integrity, limited horizontal scaling (acceptable for current scale)

## ADR-002: Adopt Celery for Async Task Processing
**Date:** 2025-12-03
**Status:** Accepted
**Context:** Batch processing requires background job execution
**Decision:** Celery with Redis broker for distributed task queue
**Consequences:** Battle-tested solution, adds Redis dependency, strong community support
```

---

### DEPENDENCIES.md

```markdown
# External Dependencies

## Third-Party Services
- **Payment Processing:** Stripe API v2023-10-16
- **Email Delivery:** SendGrid API v3
- **Monitoring:** Datadog agent for metrics and APM
- **Authentication:** Auth0 for SSO (enterprise customers only)

## Major Libraries
- **FastAPI 0.104.1:** Web framework for REST API
- **Celery 5.3.4:** Distributed task queue
- **SQLAlchemy 2.0.23:** ORM and database toolkit
- **Pydantic 2.5.0:** Data validation using Python type hints
```

---

### SECURITY.md

```markdown
# Security Architecture

## Authentication & Authorization
- **User Authentication:** JWT tokens with 1-hour expiration, refresh tokens for 30 days
- **Service-to-Service:** API keys with IP allowlisting for internal services
- **Multi-Factor Auth:** TOTP-based MFA required for admin accounts

## Data Protection
- **At Rest:** AES-256 encryption for S3 objects, encrypted RDS volumes
- **In Transit:** TLS 1.3 enforced for all external connections
- **PII Handling:** Email addresses hashed before logging, GDPR compliance

## Threat Model
- **Primary Threats:** Unauthorized data access, batch injection attacks, DDoS
- **Mitigations:** Rate limiting (100 req/min), input validation, WAF rules
```

---

### TESTING.md

```markdown
# Testing Strategy

## Test Pyramid
- **Unit Tests:** 450+ tests, 85% code coverage, run in <30s
- **Integration Tests:** 120 tests covering API endpoints and database interactions
- **End-to-End Tests:** 25 critical user journey tests in staging environment

## Testing Architecture
- **Framework:** pytest for all Python tests
- **Test Database:** PostgreSQL in Docker container, fresh schema per test class
- **Mocking Strategy:** Mock external services only (Stripe, SendGrid), use real implementations for internal components
- **CI/CD:** Tests run on every PR, must pass before merge, coverage must not decrease
```

---

## Usage Guidelines

**When to update each file:**

- Make updates immediately after implementing architectural changes
- Keep documentation in sync with code through PR process
- Review and update quarterly even if no major changes
- Use `/update-architecture` skill to ensure consistency and completeness

**Documentation principles:**

- Write for new team members joining the project
- Include diagrams for complex flows (use Mermaid)
- Keep historical context in DECISIONS.md (never delete ADRs)
- Link to code files with absolute paths when referencing implementations
