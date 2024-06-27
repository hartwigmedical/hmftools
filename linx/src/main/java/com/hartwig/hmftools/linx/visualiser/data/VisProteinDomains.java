package com.hartwig.hmftools.linx.visualiser.data;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.visualiser.file.VisProteinDomain.PROTEIN_DOMAIN_UTR;

import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.linx.visualiser.file.VisFusion;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;
import com.hartwig.hmftools.linx.visualiser.file.VisProteinDomain;

public class VisProteinDomains
{
    public static List<VisProteinDomain> readProteinDomains(final String fileName, final List<VisFusion> fusions) throws IOException
    {
        List<VisProteinDomain> all = VisProteinDomain.read(fileName);
        return proteinDomainsInFusionGenes(fusions, all);
    }

    public static List<VisProteinDomain> exonicProteinDomains(final List<VisProteinDomain> proteinDomains, final List<VisGeneExon> exons)
    {
        final List<VisProteinDomain> result = Lists.newArrayList();
        for(VisProteinDomain proteinDomain : proteinDomains)
        {
            if(proteinDomain.name().equals(PROTEIN_DOMAIN_UTR))
            {
                final List<VisGeneExon> proteinDomainExons = exons.stream()
                        .filter(x -> x.Transcript.equals(proteinDomain.Transcript)).sorted().collect(Collectors.toList());

                for(VisGeneExon exon : proteinDomainExons)
                {
                    if(proteinDomain.overlaps(exon))
                    {
                        VisProteinDomain exonicProteinDomain = VisProteinDomain.from(proteinDomain);
                        exonicProteinDomain.Start = max(exonicProteinDomain.Start, exon.start());
                        exonicProteinDomain.End = min(exonicProteinDomain.End, exon.end());
                        result.add(exonicProteinDomain);
                    }
                }
            }
            else
            {
                result.add(proteinDomain);
            }

        }

        return result;
    }

    private static List<VisProteinDomain> proteinDomainsInFusionGenes(final List<VisFusion> fusions, final List<VisProteinDomain> proteinDomains)
    {
        final Set<String> transcripts = Sets.newHashSet();
        fusions.forEach(x ->
        {
            transcripts.add(x.TranscriptUp);
            transcripts.add(x.TranscriptDown);
        });

        return proteinDomains.stream().filter(x -> transcripts.contains(x.Transcript)).collect(Collectors.toList());
    }
}
