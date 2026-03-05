# FINDING Datamodel

The FINDING datamodel module contains all the datamodels for downstream tools to process pipeline results.

The goal of this module is to make it easy to parse the datamodel:

- The module has as few dependencies as possible, to make it easy for downstream tools to depend on it.
- The module should be stable and present a clean datamodel to depend on.

## Record and RecordBuilder

The finding datamodel uses java `record`s to represent the data. [RecordBuilder](https://github.com/Randgalt/record-builder) is used to
generate builder for each record.

## Nullable and NotNull

Finding datamodel uses `jakarta.validation` for `NotNull` and `jspecify` for `Nullable`. This strange combination of annotations is
required because downstream springboot projects are only able to use `jakarta.validation`. This is a temporary solution until jspecify
is adopted by springboot.
