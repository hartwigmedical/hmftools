package com.hartwig.hmftools.linx.visualiser.circos;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.linx.visualiser.data.Exon;
import com.hartwig.hmftools.linx.visualiser.data.FusedExon;
import com.hartwig.hmftools.linx.visualiser.data.FusedExons;
import com.hartwig.hmftools.linx.visualiser.data.Fusion;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomains;

import org.jetbrains.annotations.NotNull;

public class FusionsDataWriter
{
    private final String filePrefix;

    public FusionsDataWriter(@NotNull final String sample, @NotNull final String outputDir)
    {
        this.filePrefix = outputDir + File.separator + sample;
    }

    public void write(@NotNull final List<Fusion> fusions, @NotNull final List<Exon> exons,
            @NotNull final List<ProteinDomain> proteinDomains) throws IOException
    {
        final List<FusedExon> finalExons = Lists.newArrayList();
        final List<ProteinDomain> finalProteinDomains = Lists.newArrayList();

        for (Fusion fusion : fusions)
        {
            final List<FusedExon> fusedExons = FusedExons.fusedExons(fusion, exons);
            final List<ProteinDomain> fusedProteinDomain = ProteinDomains.fusedProteinDomains(fusion, fusedExons, proteinDomains);

            final List<GenomePosition> unadjustedPositions = Lists.newArrayList();
            for (FusedExon fusedExon : fusedExons)
            {
                unadjustedPositions.add(GenomePositions.create(fusedExon.chromosome(), fusedExon.geneStart()));
                unadjustedPositions.add(GenomePositions.create(fusedExon.chromosome(), fusedExon.geneEnd()));
                unadjustedPositions.add(GenomePositions.create(fusedExon.chromosome(), fusedExon.start()));
                unadjustedPositions.add(GenomePositions.create(fusedExon.chromosome(), fusedExon.end()));
            }

            for (ProteinDomain proteinDomain : fusedProteinDomain)
            {
                unadjustedPositions.add(GenomePositions.create(proteinDomain.chromosome(), proteinDomain.start()));
                unadjustedPositions.add(GenomePositions.create(proteinDomain.chromosome(), proteinDomain.end()));
            }

            final ScalePosition scalePosition = new ScalePosition(unadjustedPositions);
            finalExons.addAll(scalePosition.scaleFusedExon(fusedExons));
            finalProteinDomains.addAll(scalePosition.interpolateProteinDomains(fusedProteinDomain));

        }

        FusedExons.write(filePrefix + ".fusions.tsv", finalExons);

        final ProteinDomainColors domainColors =
                new ProteinDomainColors(proteinDomains.stream().map(ProteinDomain::name).collect(Collectors.toSet()));
        ProteinDomains.write(filePrefix + ".protein_domains.tsv", domainColors, finalProteinDomains);
    }

}
