package com.hartwig.hmftools.datamodel.rna;

import org.immutables.value.Value;

@Value.Immutable
public interface GeneExpression {
    String geneName();
    double tpm();
    double medianTpmCancer();
    double percentileCancer();
    double medianTpmCohort();
    double percentileCohort();
}
