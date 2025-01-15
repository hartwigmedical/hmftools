####
# SQL updates for Pipeline release 6.1
# NOTE: only add updates to this script if the tools impacted by them will be released with this pipeline release

# Remove Gridss SV columns

ALTER TABLE structuralVariant
    DROP COLUMN imprecise,
    DROP COLUMN recovered,
    DROP COLUMN recoveryMethod,
    DROP COLUMN recoveryFilter,
    DROP COLUMN startRefContext,
    DROP COLUMN endRefContext;

