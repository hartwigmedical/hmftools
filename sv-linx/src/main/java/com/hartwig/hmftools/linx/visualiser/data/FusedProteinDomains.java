package com.hartwig.hmftools.linx.visualiser.data;

import static com.hartwig.hmftools.linx.visualiser.data.FusedExons.convertRegion;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.linx.visualiser.circos.ColorPicker;
import com.hartwig.hmftools.linx.visualiser.circos.ProteinDomainColors;

import org.jetbrains.annotations.NotNull;

public class FusedProteinDomains
{

    private static final String DELIMITER = "\t";

    @NotNull
    public static List<ProteinDomain> fusedProteinDomains(@NotNull final Fusion fusion, @NotNull final List<FusedExon> fusedExons,
            @NotNull final List<ProteinDomain> proteinDomains)
    {
        final List<ProteinDomain> result = Lists.newArrayList();
        if (fusedExons.isEmpty())
        {
            return result;
        }

        final FusedExon firstUpExon = fusedExons.get(0);
        final GenomeRegion upGeneRegion = upGeneRegion(fusion, firstUpExon);

        final FusedExon finalDownExon = fusedExons.get(fusedExons.size() - 1);
        final GenomeRegion downGeneRegion = downGeneRegion(fusion, finalDownExon);

        for (ProteinDomain unadjustedDomain : proteinDomains)
        {
            if (unadjustedDomain.transcript().equals(fusion.transcriptUp()) && unadjustedDomain.overlaps(upGeneRegion))
            {
                final GenomeRegion convertedDomain = convertRegion(fusion.strandUp(), upGeneRegion, unadjustedDomain);
                final ProteinDomain domain = ImmutableProteinDomain.builder().from(unadjustedDomain)
                        .chromosome(fusion.name())
                        .start(Math.max(convertedDomain.start(), firstUpExon.geneStart()))
                        .end(Math.min(convertedDomain.end(), firstUpExon.geneEnd()))
                        .build();

                result.add(domain);
            }

            if (unadjustedDomain.transcript().equals(fusion.transcriptDown()) && unadjustedDomain.overlaps(downGeneRegion))
            {
                final GenomeRegion convertedDomain = convertRegion(fusion.strandDown(), downGeneRegion, unadjustedDomain);
                final ProteinDomain domain = ImmutableProteinDomain.builder().from(unadjustedDomain)
                        .chromosome(fusion.name())
                        .start(Math.max(convertedDomain.start() + firstUpExon.geneEnd(), finalDownExon.geneStart()))
                        .end(Math.min(convertedDomain.end() + firstUpExon.geneEnd(), finalDownExon.geneEnd()))
                        .build();

                result.add(domain);
            }
        }

        return result;
    }

    @NotNull
    private static GenomeRegion upGeneRegion(@NotNull final Fusion fusion, @NotNull final FusedExon firstUpGene)
    {
        final long upGeneLength = firstUpGene.geneEnd() - firstUpGene.geneStart();
        final long upGeneStart = fusion.strandUp() < 0 ? fusion.positionUp() : fusion.positionUp() - upGeneLength;

        return GenomeRegions.create(fusion.chromosomeUp(), upGeneStart, upGeneStart + upGeneLength);
    }

    @NotNull
    private static GenomeRegion downGeneRegion(@NotNull final Fusion fusion, @NotNull final FusedExon finalDownExon)
    {
        final long downGeneLength = finalDownExon.geneEnd() - finalDownExon.geneStart();
        final long downGeneStart = fusion.strandDown() < 0 ? fusion.positionDown() - downGeneLength : fusion.positionDown();

        return GenomeRegions.create(fusion.chromosomeDown(), downGeneStart, downGeneStart + downGeneLength);
    }

    public static void write(@NotNull final String fileName, @NotNull final ProteinDomainColors colors,
            @NotNull final List<ProteinDomain> domains) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(colors, domains));
    }

    @NotNull
    static List<String> toLines(@NotNull final ProteinDomainColors colors, @NotNull final List<ProteinDomain> domains)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        domains.stream().map(x -> toString(colors, x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER).add("sampleId")
                .add("clusterId")
                .add("fusion")
                .add("start")
                .add("end")
                .add("name")
                .add("color")
                .add("transcript")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final ProteinDomainColors colors, @NotNull final ProteinDomain domain)
    {
        return new StringJoiner(DELIMITER)
                .add(domain.sampleId())
                .add(String.valueOf(domain.clusterId()))
                .add(domain.chromosome())
                .add(String.valueOf(domain.start()))
                .add(String.valueOf(domain.end()))
                .add(domain.name())
                .add(ColorPicker.hexColor(colors.color(domain)))
                .add(domain.transcript())
                .toString();
    }

}
