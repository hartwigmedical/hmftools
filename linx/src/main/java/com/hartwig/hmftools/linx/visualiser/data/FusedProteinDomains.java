package com.hartwig.hmftools.linx.visualiser.data;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
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
import com.hartwig.hmftools.linx.visualiser.file.VisFusion;
import com.hartwig.hmftools.linx.visualiser.file.VisProteinDomain;

import org.jetbrains.annotations.NotNull;

public class FusedProteinDomains
{
    public static List<VisProteinDomain> fusedProteinDomains(final VisFusion fusion, final List<FusedExon> fusedExons,
            final List<VisProteinDomain> proteinDomains)
    {
        final List<VisProteinDomain> result = Lists.newArrayList();
        if(fusedExons.isEmpty())
        {
            return result;
        }

        final FusedExon firstUpExon = fusedExons.get(0);
        final GenomeRegion upGeneRegion = upGeneRegion(fusion, firstUpExon);

        final FusedExon finalDownExon = fusedExons.get(fusedExons.size() - 1);
        final GenomeRegion downGeneRegion = downGeneRegion(fusion, finalDownExon);

        for(VisProteinDomain unadjustedDomain : proteinDomains)
        {
            if(unadjustedDomain.Transcript.equals(fusion.TranscriptUp) && unadjustedDomain.overlaps(upGeneRegion))
            {
                final GenomeRegion convertedDomain = convertRegion(fusion.StrandUp, upGeneRegion, unadjustedDomain);

                final VisProteinDomain domain = VisProteinDomain.from(unadjustedDomain);
                domain.Chromosome = fusion.name();
                domain.Start = max(convertedDomain.start(), firstUpExon.geneStart());
                domain.End = min(convertedDomain.end(), firstUpExon.geneEnd());
                result.add(domain);
            }

            if(unadjustedDomain.Transcript.equals(fusion.TranscriptDown) && unadjustedDomain.overlaps(downGeneRegion))
            {
                final GenomeRegion convertedDomain = convertRegion(fusion.StrandDown, downGeneRegion, unadjustedDomain);

                final VisProteinDomain domain = VisProteinDomain.from(unadjustedDomain);
                domain.Chromosome = fusion.name();
                domain.Start = max(convertedDomain.start() + firstUpExon.geneEnd(), finalDownExon.geneStart());
                domain.End = min(convertedDomain.end() + firstUpExon.geneEnd(), finalDownExon.geneEnd());
                result.add(domain);
            }
        }

        return result;
    }

    private static GenomeRegion upGeneRegion(final VisFusion fusion, final FusedExon firstUpGene)
    {
        final int upGeneLength = firstUpGene.geneEnd() - firstUpGene.geneStart();
        final int upGeneStart = fusion.StrandUp < 0 ? fusion.PosUp : fusion.PosUp - upGeneLength;

        return GenomeRegions.create(fusion.ChrUp, upGeneStart, upGeneStart + upGeneLength);
    }

    private static GenomeRegion downGeneRegion(final VisFusion fusion, final FusedExon finalDownExon)
    {
        final int downGeneLength = finalDownExon.geneEnd() - finalDownExon.geneStart();
        final int downGeneStart = fusion.StrandDown < 0 ? fusion.PosDown - downGeneLength : fusion.PosDown;

        return GenomeRegions.create(fusion.ChrDown, downGeneStart, downGeneStart + downGeneLength);
    }

    public static void write(final String fileName, final ProteinDomainColors colors,
            final List<VisProteinDomain> domains) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(colors, domains));
    }

    static List<String> toLines(final ProteinDomainColors colors, final List<VisProteinDomain> domains)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        domains.stream().map(x -> toString(colors, x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add("sampleId")
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
    private static String toString(final ProteinDomainColors colors, final VisProteinDomain domain)
    {
        return new StringJoiner(TSV_DELIM)
                .add(domain.SampleId)
                .add(String.valueOf(domain.ClusterId))
                .add(domain.chromosome())
                .add(String.valueOf(domain.start()))
                .add(String.valueOf(domain.end()))
                .add(domain.name())
                .add(ColorPicker.hexColor(colors.color(domain)))
                .add(domain.Transcript)
                .toString();
    }

}
