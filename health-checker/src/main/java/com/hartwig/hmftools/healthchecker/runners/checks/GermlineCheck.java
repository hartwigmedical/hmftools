package com.hartwig.hmftools.healthchecker.runners.checks;

public enum GermlineCheck {
    GERMLINE_INDEL_COUNT,
    GERMLINE_SNP_COUNT,
    GERMLINE_HETEROZYGOUS_COUNT,
    GERMLINE_HETEROZYGOUS_COUNT_ABOVE_50_VAF,
    GERMLINE_HETEROZYGOUS_COUNT_BELOW_50_VAF
}
