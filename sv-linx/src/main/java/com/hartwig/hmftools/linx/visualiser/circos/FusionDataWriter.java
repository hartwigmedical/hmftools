package com.hartwig.hmftools.linx.visualiser.circos;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.linx.visualiser.data.Exon;
import com.hartwig.hmftools.linx.visualiser.data.FusedExon;
import com.hartwig.hmftools.linx.visualiser.data.FusedExons;
import com.hartwig.hmftools.linx.visualiser.data.FusedProteinDomains;
import com.hartwig.hmftools.linx.visualiser.data.Fusion;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.VisProteinDomains;

import org.jetbrains.annotations.NotNull;

public class FusionDataWriter
{
    private final List<FusedExon> finalExons;
    private final List<ProteinDomain> finalProteinDomains;
    private final ProteinDomainColors proteinDomainColors;

    public FusionDataWriter(@NotNull final List<Fusion> fusions, @NotNull final List<Exon> exons,
            @NotNull final List<ProteinDomain> proteinDomains)
    {

        this.finalExons = Lists.newArrayList();
        this.finalProteinDomains = Lists.newArrayList();

        final List<ProteinDomain> exonicProteinDomains = VisProteinDomains.exonicProteinDomains(proteinDomains, exons);

        for (Fusion fusion : fusions)
        {
            final List<FusedExon> fusedExons = FusedExons.fusedExons(fusion, exons);
            final List<ProteinDomain> fusedProteinDomain = FusedProteinDomains.fusedProteinDomains(fusion, fusedExons, exonicProteinDomains);

            final ScaleIntrons scaler = new ScaleIntrons(ScaleIntrons.introns(fusedExons));
            final List<FusedExon> intronScaledExons = scaler.scaleIntronsInExons(fusedExons);
            final List<ProteinDomain> intronScaledProteinDomains = scaler.scaleIntronsInProteinDomains(fusedProteinDomain);

            final List<GenomePosition> unadjustedPositions = Lists.newArrayList();
            for (FusedExon fusedExon : intronScaledExons)
            {
                unadjustedPositions.add(GenomePositions.create(fusedExon.fusion(), fusedExon.geneStart()));
                unadjustedPositions.add(GenomePositions.create(fusedExon.fusion(), fusedExon.geneEnd()));
                unadjustedPositions.add(GenomePositions.create(fusedExon.fusion(), fusedExon.start()));
                unadjustedPositions.add(GenomePositions.create(fusedExon.fusion(), fusedExon.end()));
            }

            for (ProteinDomain proteinDomain : intronScaledProteinDomains)
            {
                unadjustedPositions.add(GenomePositions.create(proteinDomain.chromosome(), proteinDomain.start()));
                unadjustedPositions.add(GenomePositions.create(proteinDomain.chromosome(), proteinDomain.end()));
            }

            final ScalePosition scalePosition = new ScalePosition(unadjustedPositions);
            finalExons.addAll(scalePosition.scaleFusedExon(intronScaledExons));
            finalProteinDomains.addAll(scalePosition.interpolateProteinDomains(intronScaledProteinDomains));
        }

        this.proteinDomainColors = new ProteinDomainColors(finalProteinDomains);
    }

    @NotNull
    public Object write(@NotNull final String sample, @NotNull final String outputDir)
            throws IOException
    {
        String filePrefix = outputDir + File.separator + sample;
        FusedExons.write(filePrefix + ".fusions.tsv", finalExons);
        FusedProteinDomains.write(filePrefix + ".protein_domains.tsv", proteinDomainColors, finalProteinDomains);
        return this;
    }

    public List<FusedExon> finalExons()
    {
        return finalExons;
    }

}
