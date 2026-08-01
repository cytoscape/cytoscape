# Contributing to Cytoscape

Cytoscape is an open-source platform for network analysis and visualization. We appreciate contributions of all kinds, including code, tests, documentation, bug reports, feature proposals, and usability improvements.

Cytoscape Desktop is composed of multiple repositories. The [`cytoscape/cytoscape`](https://github.com/cytoscape/cytoscape) repository is the top-level Maven project used to assemble the complete Cytoscape Desktop distribution.

Before beginning a substantial change, please open an issue or discuss the proposal with the Cytoscape development team. This helps confirm that the change belongs in the correct repository and is consistent with the project’s architecture and roadmap.

## Ways to contribute

You can contribute by:

- Reporting reproducible bugs
- Fixing existing issues
- Adding or improving automated tests
- Improving documentation
- Proposing or implementing features
- Improving accessibility, usability, or performance
- Reviewing pull requests
- Developing or maintaining Cytoscape Apps

Changes to independently released Cytoscape Apps should normally be submitted to the App’s own repository rather than this top-level repository.

## Reporting issues

Before opening an issue:

1. Search the existing issues to check whether the problem has already been reported.
2. Confirm that the issue occurs in a currently supported Cytoscape version or recent development build.
3. Determine which Cytoscape repository contains the affected component, when possible.

A useful bug report includes:

- Cytoscape version or commit
- Operating system and version
- Java version
- Steps needed to reproduce the problem
- Expected behavior
- Actual behavior
- Relevant logs, screenshots, sample networks, or session files
- Whether the problem occurs with third-party Apps disabled

Do not include confidential or sensitive data in public issues or attachments.

Feature requests should explain the use case and the problem being solved, rather than describing only a proposed implementation.

## Choosing the correct repository

Cytoscape Desktop is a multi-repository project. The main components include:

- [`cytoscape/cytoscape`](https://github.com/cytoscape/cytoscape): Top-level build and workspace-management project
- [`cytoscape/cytoscape-api`](https://github.com/cytoscape/cytoscape-api): Public Cytoscape APIs
- [`cytoscape/cytoscape-impl`](https://github.com/cytoscape/cytoscape-impl): Core implementations
- [`cytoscape/cytoscape-support`](https://github.com/cytoscape/cytoscape-support): Supporting bundles and utilities
- [`cytoscape/cytoscape-gui-distribution`](https://github.com/cytoscape/cytoscape-gui-distribution): Desktop distribution and assembly
- Individual core-App repositories in the [`cytoscape`](https://github.com/cytoscape) GitHub organization

Open the issue and pull request in the repository that owns the code being changed. When a change spans several repositories, describe the dependencies between the corresponding pull requests and link them to one another.

## Development requirements

Building the current Cytoscape Desktop development version requires:

- Git
- JDK 17
- Maven 3
- A Unix-compatible shell for `cy.sh`, or an equivalent environment on Windows

Confirm the installed versions with:

```sh
git --version
java -version
mvn --version
```

Set `JAVA_HOME` to the JDK 17 installation when necessary.

## Setting up the source tree

Fork [`cytoscape/cytoscape`](https://github.com/cytoscape/cytoscape), then clone your fork:

```sh
git clone https://github.com/YOUR-USERNAME/cytoscape.git
cd cytoscape
```

Add the upstream repository:

```sh
git remote add upstream https://github.com/cytoscape/cytoscape.git
git fetch upstream
```

Clone the Cytoscape core repositories. Run this from inside your clone of this
repository — the script lives here, and it clones the core subprojects into this
same directory, alongside `pom.xml`:

```sh
./cy.sh pull
```

Re-run it any time to update every repository; it clones what is missing and
pulls what is already there.

The README has a table of the `cy.sh` commands — each one has a core-repository
form and a core-app form.

## Branches

Core Cytoscape development takes place on the `develop` branch. The `master` branch is used for released code.

Create your working branch from the latest upstream `develop` branch:

```sh
git fetch upstream
git switch develop
git reset --hard upstream/develop
git switch -c issue-123-short-description
```

Use a focused branch for each change. Avoid combining unrelated fixes or refactoring with a feature contribution.

Core Apps may use a different branching model. Check the target App repository before creating a branch.

## Building Cytoscape

After initializing the workspace, enter the generated project directory:

```sh
cd cytoscape
```

Build the complete project with:

```sh
mvn clean install -U
```

For a first build, do not skip tests. Some modules depend on test outputs produced by earlier modules.

After a successful build, the assembled application is located under:

```text
gui-distribution/assembly/target/cytoscape
```

Run it on macOS or Linux with:

```sh
./gui-distribution/assembly/target/cytoscape/cytoscape.sh
```

Run it on Windows with:

```bat
gui-distribution\assembly\target\cytoscape\cytoscape.bat
```

For information about platform-specific setup or known Maven issues, consult the main repository’s [`README.md`](https://github.com/cytoscape/cytoscape/blob/develop/README.md).

## Making changes

Follow the existing architecture and conventions of the module you are modifying.

In particular:

- Keep public interfaces in the appropriate API modules.
- Keep implementations in implementation modules.
- Preserve compatibility unless an API change has been explicitly discussed and approved.
- Avoid introducing dependencies between modules without a clear architectural need.
- Keep commits and pull requests narrowly scoped.
- Do not commit generated build output, IDE metadata, local configuration, or credentials.
- Update documentation when behavior, APIs, configuration, or user-visible functionality changes.

Public API changes require additional review because they can affect Cytoscape Apps and other downstream integrations. Discuss proposed API additions or incompatible changes before implementation.

## Java code style

Match the style of the surrounding code.

As general guidance:

- Use tabs or spaces consistently with the existing module.
- Use descriptive class, method, and variable names.
- Keep methods focused and avoid unnecessary complexity.
- Prefer interfaces and services over tightly coupled implementations.
- Add Javadoc to public APIs and non-obvious behavior.
- Remove unused imports and commented-out code.
- Avoid unrelated formatting changes.
- Handle errors explicitly rather than silently ignoring them.
- Follow the existing OSGi service and bundle conventions.

A contribution should optimize for maintainability and clarity rather than minimizing the number of lines changed.

## Testing

Add or update tests for behavior affected by your change.

For a bug fix, include a regression test that fails without the fix whenever practical. For a new feature, test normal behavior, boundary cases, and relevant failure conditions.

Run the tests for the changed module:

```sh
mvn clean test
```

Then build the module and its required dependencies:

```sh
mvn clean install
```

Before submitting a pull request, build the complete Cytoscape distribution from the workspace root:

```sh
mvn clean install -U
```

Some desktop behavior cannot be adequately covered by unit tests. For user-interface changes, also test the assembled application manually and describe the scenarios tested in the pull request.

Consider testing with a clean Cytoscape configuration when existing settings or installed Apps could affect the result. Back up any configuration you need before removing or renaming it.

## Core Apps

Core Apps are maintained in separate repositories. To clone them all into an
`apps` directory inside your clone of this project — and to update them later —
run:

```sh
./cy.sh pull-apps
```

Each command that acts on the core repositories has an apps counterpart:
`pull-apps`, `switch-apps`, `branch-apps` and `build-apps`. Note that Core Apps
use `master` rather than `develop`.

An individual App may also be built on its own from its directory:

```sh
mvn clean install -U
```

A locally built App JAR can be tested by:

- Selecting **Apps → App Store → Install Apps from File** in Cytoscape, or
- Copying the JAR to the appropriate Cytoscape configuration directory

A change to a core App should generally be submitted to that App’s repository. A separate change to `cytoscape-gui-distribution` may be required when updating the version included in the desktop distribution.

## Documentation

Update documentation as part of any change that affects:

- Public APIs
- User-visible behavior
- Build or configuration procedures
- Commands or Automation interfaces
- App-development guidance
- Compatibility requirements

User-manual changes may belong in the [`cytoscape/cytoscape-manual`](https://github.com/cytoscape/cytoscape-manual) repository rather than the source repository.

Documentation should describe the current behavior directly and should not rely solely on issue or pull-request discussions.

## Commit messages

Write clear commit messages that explain the purpose of the change.

A useful subject line:

- Uses the imperative mood
- Is concise
- Identifies the affected component when helpful
- Does not end with a period

For example:

```text
Fix network view disposal in presentation manager
```

Use the commit body to explain why the change is needed, important implementation decisions, and any compatibility considerations.

Reference the related issue when applicable:

```text
Fixes #123
```

## Submitting a pull request

Push your branch to your fork and open a pull request against the appropriate Cytoscape repository.

For Cytoscape core repositories, target the `develop` branch unless a maintainer has requested another branch.

A pull request should include:

- A concise description of the problem
- A summary of the solution
- A link to the related issue
- Tests added or updated
- Manual testing performed
- Screenshots or recordings for user-interface changes
- Documentation changes
- Compatibility or migration considerations
- Links to related pull requests in other Cytoscape repositories

Keep the pull request focused. Large changes are easier to review when divided into independently meaningful commits or a coordinated series of pull requests.

By submitting a pull request, you agree that your contribution may be distributed under the repository’s license.

## Review process

Maintainers may request changes related to correctness, tests, architecture, compatibility, documentation, or maintainability.

During review:

- Respond to comments or mark resolved discussions only after addressing them.
- Push revisions to the existing pull-request branch rather than opening a replacement pull request.
- Avoid force-pushing after review has begun unless necessary.
- Rebase or merge the latest target branch when requested.
- Keep automated checks passing.

Approval does not guarantee immediate merge. Changes involving public APIs, release coordination, or several repositories may require additional review.

## Security issues

Do not publicly disclose a vulnerability before maintainers have had an opportunity to investigate and prepare a fix.

Report suspected security vulnerabilities privately through GitHub’s security-reporting mechanism when it is enabled for the affected repository, or contact the Cytoscape team through an appropriate private project channel.

## Community standards

Be respectful and constructive in issues, pull requests, reviews, and project discussions.

Critique code and technical decisions rather than individuals. Contributors are expected to follow the project’s Code of Conduct wherever one is provided.

Thank you for helping improve Cytoscape.
