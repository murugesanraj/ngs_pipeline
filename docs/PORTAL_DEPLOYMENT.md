# Secure portal deployment

## Recommended model: one Open OnDemand session per user

Hellbender publicly provides Open OnDemand and RStudio Server. In the
recommended model, each client:

1. authenticates to Hellbender with their own University account;
2. starts an Open OnDemand/RStudio compute allocation;
3. launches the Shiny app as that Unix user;
4. sees only manifests/results readable through approved group membership or
   ACLs.

This keeps authentication and authorization at institutional boundaries. The
app has no password table, role editor, or impersonation feature.

## Publishing workflow

1. Complete technical QC and scientific review.
2. Confirm that sample IDs and portal description contain no prohibited client
   or participant identifiers.
3. Confirm the correct project group/ACL and remove unintended broad read bits.
4. Run `bin/publish_project.sh PROJECT_DIR PORTAL_ROOT`.
5. Test using a non-admin client account from the intended group.
6. Record the release/version in the client handoff.

The publish command creates only a symlink. It never calls recursive `chmod` or
`chgrp`, because access changes require an explicit project-level decision.

## App controls

- Project discovery is one directory level under `NGS_PORTAL_ROOT`.
- A project must have a readable pipeline-generated manifest.
- Files must be explicitly listed in that manifest.
- Canonical path checks prevent `..`, absolute paths, and symlink escape.
- Downloads are read-only and size-limited (250 MB by default).
- MultiQC HTML is downloaded, not embedded, to avoid executing project HTML in
  the portal page.
- Raw FASTQ, BAM, and private configuration paths are not included by default.

## Central service model

A persistent, multi-user URL is an infrastructure service, not just an R
script. Before central deployment, coordinate with ITRSS for:

- approved RNet host and reverse proxy;
- University SSO and per-request identity propagation;
- TLS certificates;
- per-user filesystem authorization or a reviewed access-control service;
- audit logging, patching, uptime, monitoring, and incident response;
- accessibility review and data-classification approval.

Do not bind the app to `0.0.0.0` on a compute node or treat a shared password as
institutional authentication.

## Client access test

For every release, verify:

- intended client can open the project and approved small downloads;
- unintended HPC user cannot read the manifest or target result files;
- client cannot navigate to another project/path;
- no identifying FASTQ paths or private config appear in the UI;
- stale/withdrawn projects are removed from the portal root;
- project data remain available from the source allocation after the symlink is
  created.

