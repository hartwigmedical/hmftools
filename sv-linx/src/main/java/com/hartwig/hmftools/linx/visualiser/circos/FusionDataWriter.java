package com.hartwig.hmftools.linx.visualiser.circos;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.linx.visualiser.data.Exon;
import com.hartwig.hmftools.linx.visualiser.data.FusedExon;
import com.hartwig.hmftools.linx.visualiser.data.FusedExons;
import com.hartwig.hmftools.linx.visualiser.data.FusedProteinDomains;
import com.hartwig.hmftools.linx.visualiser.data.Fusion;
import com.hartwig.hmftools.linx.visualiser.data.ImmutableProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;

import org.jetbrains.annotations.NotNull;

public class FusionDataWriter
{
    private final List<FusedExon> finalExons;
    private final List<ProteinDomain> finalProteinDomains;
    private final ProteinDomainColors proteinDomainColors;

    public FusionDataWriter(@NotNull final List<Fusion> fusions,
            @NotNull final List<Exon> exons,
            @NotNull final List<ProteinDomain> proteinDomains, @NotNull final ProteinDomainColors proteinDomainColors)
    {
        this.finalExons = Lists.newArrayList();
        this.finalProteinDomains = Lists.newArrayList();
        this.proteinDomainColors = proteinDomainColors;

        for (Fusion fusion : fusions)
        {
            final List<FusedExon> fusedExons = FusedExons.fusedExons(fusion, exons);
            final List<ProteinDomain> fusedProteinDomain = FusedProteinDomains.fusedProteinDomains(fusion, fusedExons, proteinDomains);

            final List<GenomePosition> unadjustedPositions = Lists.newArrayList();
            for (FusedExon fusedExon : fusedExons)
            {
                unadjustedPositions.add(GenomePositions.create(fusedExon.fusion(), fusedExon.geneStart()));
                unadjustedPositions.add(GenomePositions.create(fusedExon.fusion(), fusedExon.geneEnd()));
                unadjustedPositions.add(GenomePositions.create(fusedExon.fusion(), fusedExon.start()));
                unadjustedPositions.add(GenomePositions.create(fusedExon.fusion(), fusedExon.end()));
            }

            for (ProteinDomain proteinDomain : fusedProteinDomain)
            {
                unadjustedPositions.add(GenomePositions.create(proteinDomain.chromosome(), proteinDomain.start()));
                unadjustedPositions.add(GenomePositions.create(proteinDomain.chromosome(), proteinDomain.end()));
            }

            final ScalePosition scalePosition = new ScalePosition(unadjustedPositions);
            finalExons.addAll(scalePosition.scaleFusedExon(fusedExons));
            finalProteinDomains.addAll(scalePosition.interpolateProteinDomains(fusedProteinDomain));
            finalProteinDomains.addAll(legendOnlyDomains(fusion.name(), proteinDomains, fusedProteinDomain));
        }
    }

    public void write(@NotNull final String sample, @NotNull final String outputDir)
            throws IOException
    {
        String filePrefix = outputDir + File.separator + sample;
        FusedExons.write(filePrefix + ".fusions.tsv", finalExons);
        FusedProteinDomains.write(filePrefix + ".protein_domains.tsv", proteinDomainColors, finalProteinDomains);
    }

    public List<FusedExon> finalExons()
    {
        return finalExons;
    }

    @NotNull
    private Set<ProteinDomain> legendOnlyDomains(@NotNull final String fusion, @NotNull final List<ProteinDomain> original,
            @NotNull final List<ProteinDomain> fused)
    {
        final Set<ProteinDomain> result = Sets.newHashSet();

        final Set<String> fusedNames = fused.stream().map(ProteinDomain::name).collect(Collectors.toSet());
        for (final ProteinDomain proteinDomain : original)
        {
            if (!fusedNames.contains(proteinDomain.name()))
            {
                result.add(ImmutableProteinDomain.builder().from(proteinDomain).chromosome(fusion).start(0).end(0).build());
            }
        }

        return result;
    }

}
