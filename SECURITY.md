# Security policy

## Supported version

Security fixes are applied to the latest release on the default branch.

## Reporting a concern

Do not open a public issue containing credentials, sample identifiers, client
information, protected data, internal hostnames, or storage paths. Report the
problem privately to the repository maintainer and follow University incident
reporting requirements when institutional data may be affected.

## Intended security boundary

The Shiny application is read-only and is designed to run as the authenticated
Hellbender user through Open OnDemand. POSIX group permissions/ACLs determine
which project manifests and result files that process can read. The repository
does not implement a password database or replace institutional authentication.

Do not expose the app's TCP port directly to an untrusted network. A shared,
central deployment requires an institutionally approved reverse proxy and SSO
configuration reviewed by IT Research Support Solutions.

## Data limitations

No production data, re-identification keys, secrets, or client-specific config
belongs in Git. Follow the current Hellbender/RDE data-classification policy and
the project's approved data-management plan.

