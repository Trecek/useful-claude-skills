# Investigation: Batch Processing Drops Unicode Records

> **Example Note:** This is a sample output for demonstration purposes.
> Actual outputs will reflect your specific codebase and issues.

## Summary

Investigation into reported issue where batch processing silently drops records containing Unicode characters. Root cause identified in the CSV parser's encoding assumptions. The parser defaults to ASCII encoding when reading input files, causing decode errors that are caught and suppressed by an overly broad exception handler. Records with Unicode characters fail to decode, trigger the exception, and are silently skipped. This affects approximately 3-8% of production batches based on log analysis.

## Root Cause

**File:** `src/batch/csv_parser.py`, lines 45-62

The `CSVParser.read_batch()` method opens files without specifying encoding, defaulting to the system locale (typically ASCII on production servers). When Unicode characters are encountered:

1. The file reader raises `UnicodeDecodeError`
2. The exception is caught by a broad `except Exception` block (line 58)
3. The error is logged at DEBUG level only
4. Processing continues with the next record, silently dropping the failed one

```python
def read_batch(self, filepath):
    records = []
    try:
        with open(filepath) as f:  # No encoding specified - defaults to ASCII
            for line in f:
                records.append(self._parse_line(line))
    except Exception as e:  # Too broad - catches UnicodeDecodeError
        logger.debug(f"Error processing record: {e}")  # Silent failure
    return records
```

## Affected Components

**Primary:**
- `src/batch/csv_parser.py` - CSV parsing logic with encoding issue
- `src/batch/processor.py` - Calls the parser, expects all records returned
- `src/batch/validator.py` - Validates batch completeness but relies on parsed count

**Secondary:**
- `src/api/batch_upload.py` - Accepts file uploads, could validate encoding earlier
- `src/monitoring/metrics.py` - Records processing metrics but misses dropped records
- `tests/batch/test_csv_parser.py` - Existing tests only use ASCII test data

## Data Flow

1. **Upload:** User uploads CSV via `api/batch_upload.py` endpoint
2. **Storage:** File saved to staging directory with original encoding
3. **Processing:** `batch/processor.py` invokes `csv_parser.read_batch()`
4. **Parsing:** Parser opens file (ASCII default), reads line by line
5. **Error Point:** Unicode characters cause `UnicodeDecodeError`, caught silently
6. **Validation:** `validator.py` compares parsed count against expected, but expected count is derived from the incomplete parsed data
7. **Output:** Incomplete batch written to database, appears successful

## Test Gap Analysis

**Why existing tests didn't catch this:**

1. **Test Data Homogeneity:** All test fixtures in `tests/fixtures/sample_batches/` use pure ASCII data
2. **No Edge Case Coverage:** No tests for international characters, emojis, or multi-byte encodings
3. **Mock Overuse:** Integration tests mock the file I/O layer, bypassing actual encoding issues
4. **Metrics Not Validated:** Tests verify processing completes but don't validate record counts match input
5. **Log Level in Tests:** DEBUG logs are suppressed in test runs, masking the warning signs

**Missing test categories:**
- Unicode characters (Latin extended, Cyrillic, CJK)
- Emoji and special symbols
- Mixed encoding files
- Malformed UTF-8 sequences
- Large batches with sporadic Unicode (realistic scenario)

## Similar Patterns

**Good patterns found elsewhere:**

1. **`src/reports/excel_reader.py`** (lines 23-28): Explicitly specifies `encoding='utf-8'` and validates with a BOM check
2. **`src/import/json_loader.py`** (lines 12-15): Uses encoding detection library (`chardet`) before parsing
3. **`src/export/csv_writer.py`** (lines 34-40): Documents encoding in file header comment

**Anti-patterns to avoid:**

1. **`src/legacy/data_importer.py`** (lines 89-103): Similar broad exception handling without proper logging
2. **`src/batch/xml_parser.py`** (lines 56-61): Assumes UTF-8 but doesn't handle errors gracefully

## External Research

**Unicode in CSV Processing:**
- Python's `open()` defaults to platform encoding (PEP 529, Python 3.6+)
- Best practice: Always specify `encoding='utf-8'` explicitly (Real Python, "Unicode & Character Encodings in Python", 2024)
- CSV RFC 4180 doesn't specify encoding; UTF-8 with BOM is common convention

**Error Handling Patterns:**
- Broad exception catching considered anti-pattern (PEP 8, Exception Handling Guidelines)
- Silent failures violate "errors should never pass silently" (PEP 20, Zen of Python)
- Log at WARNING level minimum for data integrity issues (Python Logging Cookbook, 2025)

**Production Impact:**
- Affected records often from international customers (EU, APAC regions)
- Similar issues reported in Apache Commons CSV (JIRA: CSV-287) and pandas (GitHub issue #15086)

## Recommendations

**Immediate fixes:**

1. **Explicit Encoding:** Change `open(filepath)` to `open(filepath, encoding='utf-8', errors='strict')`
2. **Proper Error Handling:** Replace broad `except Exception` with specific exception types
3. **Error Logging:** Raise errors to ERROR level, include record context
4. **Validation:** Add pre-processing encoding validation step

**Long-term improvements:**

1. **Encoding Detection:** Implement `chardet`-based detection for ambiguous files
2. **Early Validation:** Add encoding check at upload API before accepting files
3. **Metrics Enhancement:** Track dropped records separately, alert on anomalies
4. **Test Data Diversity:** Expand test fixtures to include realistic international data
5. **Documentation:** Add encoding requirements to API documentation and user guides

**Testing strategy:**

1. Add unit tests with Unicode test data (multiple character sets)
2. Add integration test that verifies input line count matches output record count
3. Add regression test specifically for this bug with real problematic data
4. Add property-based tests (hypothesis library) generating random Unicode strings

**Rollout considerations:**

- Backward compatibility: Existing ASCII files unaffected
- Migration: May need one-time reprocessing of recent batches
- Monitoring: Add alerts for encoding errors during transition period
