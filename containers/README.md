# Tool containers

All published tool images use one Dockerfile and one locked Pixi workspace.
Images remain independently versioned; only their shared build definition is
centralized.

- `pixi.toml` is the source of truth. Every environment is named after its Maven
  module, and its presence under `[environments]` means that module publishes an
  image.
- Pixi features hold the small dependency sets shared by related tools.
- `pixi.lock` makes every Conda/Bioconda package reproducible.
- `Dockerfile` installs one selected environment and copies only that
  environment and the selected tool artifact into the runtime image.
- Java 17 comes from the pinned official Eclipse Temurin JRE image, avoiding the
  large graphical dependency closure of Conda's OpenJDK package.

Images target `linux/amd64`, matching Bioconda package availability and the
existing release images.

Use Pixi 0.73 or newer (the version pinned in Docker and CI). To add an image,
add an environment to `pixi.toml`, reuse or add a focused feature, and refresh
the lock:

```sh
pixi lock --manifest-path containers/pixi.toml
```

CI verifies that the committed lock still matches the manifest. Container builds
use `--locked`, so releases cannot silently change dependency versions.

`health-checker/Dockerfile` was removed because that Maven module no longer
exists; QSee replaced HealthChecker.
