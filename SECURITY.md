# Security

The Flask app is a **local development interface** bound to 127.0.0.1;
do not expose it directly to a network. Uploaded coordinates may be
confidential; no cloud uploads are performed by this software itself.
The UI and analysis are not designed for untrusted public Internet users.
Inputs are limited by size; malicious structure files or browser dependencies
can still contain vulnerabilities. Report vulnerabilities privately to
the repository maintainer using GitHub's private reporting facilities when
available, not a public issue with exploit details or private files.
Never commit secrets, credentials or unpublished experimental structures.
