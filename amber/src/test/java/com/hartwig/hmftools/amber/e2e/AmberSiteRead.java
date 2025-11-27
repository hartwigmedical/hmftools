package com.hartwig.hmftools.amber.e2e;

import java.util.Set;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.bam.testutilities.BasesRegion;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

class AmberSiteRead
{
    private static final Set<String> ALL_BASES = ImmutableSet.of("A", "C", "G", "T");
    private static final int READ_SEMI_LENGTH = 5;
    private static final int READ_DISTANCE = 20;

    public enum BaseInstruction
    {
        REF
                {
                    @Override
                    public String getBase(final AmberSite site)
                    {
                        return site.ref();
                    }
                },
        ALT
                {
                    @Override
                    public String getBase(final AmberSite site)
                    {
                        return site.alt();
                    }
                },
        OTHER
                {
                    @Override
                    public String getBase(final AmberSite site)
                    {
                        Set<String> bases = Sets.newHashSet(site.ref(), site.alt());
                        return Sets.difference(ALL_BASES, bases).iterator().next();
                    }
                };

        public abstract String getBase(AmberSite site);
    }

    private final AmberSite Site;
    private final BaseInstruction Instruction;

    public AmberSiteRead(final AmberSite site, final BaseInstruction instruction)
    {
        Site = site;
        Instruction = instruction;
    }

    HumanChromosome chromosome()
    {
        return Site.chr();
    }

    public BasesRegion baseRegion(RefGenomeSource refGenomeSource)
    {
        int leftPartStart = Site.Position - (READ_SEMI_LENGTH - 1);
        byte[] leftPart = refGenomeSource.getBases(Site.Chromosome, leftPartStart, Site.Position - 1);
        int rightPartEnd = Site.Position + READ_SEMI_LENGTH;
        byte[] rightPart = refGenomeSource.getBases(Site.Chromosome, Site.Position + 1, rightPartEnd);
        byte[] bytes = new byte[2 * READ_SEMI_LENGTH];
        System.arraycopy(leftPart, 0, bytes, 0, leftPart.length);
        System.arraycopy(Instruction.getBase(Site).getBytes(), 0, bytes, READ_SEMI_LENGTH - 1, 1);
        System.arraycopy(rightPart, 0, bytes, READ_SEMI_LENGTH, rightPart.length);
        return new BasesRegion(chromosome(), leftPartStart, rightPartEnd, bytes);
    }

    public BasesRegion mateReadBasesRegion(RefGenomeSource refGenomeSource)
    {
        int start = Site.Position - (READ_SEMI_LENGTH - 1) + READ_DISTANCE;
        int end = Site.Position + READ_SEMI_LENGTH + READ_DISTANCE;
        byte[] bytes = refGenomeSource.getBases(Site.Chromosome, start, end);
        return new BasesRegion(chromosome(), start, end, bytes);
    }
}
