# Developer Notes

This document is to orient developers of the pipeline on the structure and development cycle of the pipeline.

## Repository Structure

* `config/`: Contains general configuration `yaml` files for different use cases of the pipeline which are used to populate default values in project configuration using `scripts/config_creator.py`.
* `data/`: Contains test data for unit and functional tests (folders ending in `_test_data`) and reference data for different parts of the pipeline (blacklists and TSSs in `data/references/`).
* `envs/`: Conda environment definitions for rules in the pipeline as well as the overall pipeline conda environment (`envs/atac_chip_pipeline.yaml`).
* `rules/`: Snakemake files used in `include:` statements in the top-level snakemake files.
* `scripts/`: Scripts called by snakemake rules or to help setup the pipeline.
* `tests/`: Unit tests to run with `nosetests` prior to merges and releases.
* `user_tests/`: End-to-end functional tests for users to run after downloading the pipeline and setting up the pipeline conda environment.
* The top-level snakemake files are the actual pipelines that users call:
    * `ATACseq.smk`
    * `ChIPseq_pe.smk`
    * `ChIPseq_se.smk`
    * `homer_only.smk`

## Development Cycle

### Intention

GitHub issues are a useful way to keep track of bugfixes and new features that are waiting to be implemented in the pipeline. These can be user or developer generated. It is helpful to be descriptive with issues--writing pseudo-code, making notes for potential issues, etc. Usage of the comments in an issue can help document the implementation process, or solicit comments from users or other developers.

### Branching

The two permanent branches of the repository are `master` for production code and `develop` which may be in a self-consistent state, but may not be. A third type of branch, issue branches, can be used to isolate development of bugfixes/patches and/or new features that may require multiple commits. It is helpful to number them to match the issue in GitHub issues for tracking purposes.

### Versioning

We use [semantic versioning](https://semver.org/spec/v2.0.0.html) as a way to communicate the amount of change encapsulated in releases to users.

### Implementing

#### Bugfixes and Patches

The `Z` in semantic version's `X.Y.Z`. Usually if a bug exists in the codebase, it doesn't have a test associated with it. It's probably useful, by default, to add a test confirming the bug is fixed in the appropriate test script in the `tests/` folder. Otherwise, a clear justification for not adding a test should be thought through.

#### New Features

The `Y` in semantic version's `X.Y.Z` if the feature is backwards compatible, or the `X` if the feature is not. Some considerations when adding a new feature:

* Where does the feature fit in the existing pipeline(s)?
    * What are the `input` and `output` of the rule, and how do they relate to existing rules?
    * What wildcards, are needed for the rule?
    * Is it for use in both ATAC and/or ChIP? If both, consider whether to duplicate it in the top-level snakemake files, or if it can be generalized and put in `rules/` and then `include:`d in the top-level snakemake files.
* Does it require a new conda environment for the rule(s)? Add the environment definition `yaml` to `config/` if so.
* Do the rules associated with the feature have any parameters which should be user-configurable? If so, changes to `scripts/config_creator.py` are likely needed. Consider the existing parameter style of including them verbatim in the configuration `yaml` as they would be pasted into the command.
* Do the rules require any custom-scripts?
    * If so, consider appropriate user-guardrails and associated unit tests to include in `tests/`.

### Merging

Once a bugfix/patch or new feature is implemented in its issue branch, it should be merged into `develop` after `nosetests` complete without any failures (though exceptions can be made for merges into `develop`). It is a good rule that all `nosetests` pass, without exception, prior to merging into `master` and making a release.
