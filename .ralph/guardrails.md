# Ralph Guardrails (Signs)

- Prefer small, incremental commits.
- Do not read large files wholesale; use rg/grep to locate relevant sections first.
- Add fixtures before implementing parsers; keep fixtures tiny.
- If stuck, write a short note to `.ralph/errors.log` and simplify the step.
- Never change schemas without updating tests + docs.
- Make all decisions explicit in code (no magic inference without logging).

## Core Signs

### Sign: Read Before Writing
- **Trigger**: Before modifying any file
- **Instruction**: Always read the existing file first
- **Added after**: Core principle

### Sign: Test After Changes
- **Trigger**: After any code change
- **Instruction**: Run tests to verify nothing broke
- **Added after**: Core principle

### Sign: Commit Checkpoints
- **Trigger**: Before risky changes
- **Instruction**: Commit current working state first
- **Added after**: Core principle

---

## Learned Signs

(Signs added from observed failures will appear below)

