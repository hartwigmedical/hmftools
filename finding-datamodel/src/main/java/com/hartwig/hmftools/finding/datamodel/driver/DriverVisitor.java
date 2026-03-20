package com.hartwig.hmftools.finding.datamodel.driver;

import com.hartwig.hmftools.finding.datamodel.ChromosomeArmCopyNumber;
import com.hartwig.hmftools.finding.datamodel.Disruption;
import com.hartwig.hmftools.finding.datamodel.Fusion;
import com.hartwig.hmftools.finding.datamodel.GainDeletion;
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.Virus;

public interface DriverVisitor
{
    void visit(SmallVariant smallVariant);
    void visit(GainDeletion gainDeletion);
    void visit(Fusion fusion);
    void visit(Disruption disruption);
    void visit(Virus virus);
    void visit(ChromosomeArmCopyNumber chromosomeArmCopyNumber);
}
