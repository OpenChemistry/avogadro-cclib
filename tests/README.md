# Test data

Drop a real computational chemistry output file (optionally gzip-compressed,
i.e. `something.out.gz`) into `tests/data/` and it is automatically picked up
by the corpus-wide contract tests in `test_reader.py` and `test_cli.py` --
no code changes required.

If a file cannot be parsed into a readable molecule by design (e.g. it has
no atomic coordinates at all), add its logical name (the filename with any
`.gz` suffix stripped) to `EXPECTED_UNREADABLE` in `conftest.py`, along with
a comment explaining why, and cover it with its own regression test in
`test_regressions.py` instead.
