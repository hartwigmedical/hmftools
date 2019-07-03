package com.hartwig.hmftools.linx.visualiser.data;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Ignore;

@Ignore
public class FusedExonsTest
{
    @Ignore
    public void testStuff() throws IOException
    {
        final List<Fusion> fusions = Fusions.fromFile("/Users/jon/hmf/analysis/CPCT02390005T/CPCT02390005T.linx.fusions_detailed.csv")
                .stream()
                .filter(x -> x.clusterId() == 88)
                .collect(Collectors.toList());
        final List<Exon> exons = Exons.readExons("/Users/jon/hmf/analysis/CPCT02390005T/CPCT02390005T.linx.vis_gene_exon.tsv");

        List<FusedExon> fused = FusedExons.fusedExons(fusions.get(0), exons);
        for (FusedExon fusedExon : fused)
        {
            System.out.println(fusedExon);
        }

        final List<ProteinDomain> domains = ProteinDomains.readProteinDomains("/Users/jon/hmf/analysis/CPCT02390005T/CPCT02390005T.linx.vis_protein_domain.tsv");

        List<ProteinDomain> fusedDomain = ProteinDomains.fusedProteinDomains(fusions.get(0), fused, domains);
        for (ProteinDomain proteinDomain : fusedDomain)
        {
            System.out.println(proteinDomain);
        }

    }

}
