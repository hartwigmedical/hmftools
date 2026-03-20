package com.hartwig.hmftools.finding.datamodel.driver;

import java.util.HashSet;
import java.util.Set;

import com.hartwig.hmftools.finding.datamodel.ChromosomeArmCopyNumber;
import com.hartwig.hmftools.finding.datamodel.Disruption;
import com.hartwig.hmftools.finding.datamodel.Fusion;
import com.hartwig.hmftools.finding.datamodel.GainDeletion;
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.Virus;

public class DriverGenesGetter implements DriverVisitor
{
    private final Set<String> genes = new HashSet<>();

    public static Set<String> getGenes(Driver driver)
    {
        return new DriverGenesGetter(driver).genes;
    }

    private DriverGenesGetter(Driver driver)
    {
        driver.accept(this);
    }

    @Override
    public void visit(SmallVariant smallVariant)
    {
        genes.add(smallVariant.gene());
    }

    @Override
    public void visit(GainDeletion gainDeletion)
    {
        genes.add(gainDeletion.gene());
    }

    @Override
    public void visit(Fusion fusion)
    {
        genes.add(fusion.geneUp());
        genes.add(fusion.geneDown());
    }

    @Override
    public void visit(Disruption disruption)
    {
        genes.add(disruption.gene());
    }

    @Override
    public void visit(Virus virus)
    {
    }

    @Override
    public void visit(ChromosomeArmCopyNumber chromosomeArmCopyNumber)
    {
    }
}
