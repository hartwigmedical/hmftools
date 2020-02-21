package com.hartwig.hmftools.knowledgebasegenerator;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.knowledgebasegenerator.transvar.RefVersion;
import com.hartwig.hmftools.knowledgebasegenerator.transvar.Transvar;

public class TempRunner {

    public static void main(String[] args) throws IOException, InterruptedException {
        String refFastaPath = System.getProperty("user.home") + "/hmf/refgenome/Homo_sapiens.GRCh37.GATK.illumina.fasta";
        RefVersion refVersion = RefVersion.HG19;

        Transvar transvar = new Transvar(refFastaPath, refVersion);

        List<VariantHotspot> hotspots = transvar.extractHotspotsFromProteinAnnotation("MTOR",  "L2230V");

        System.out.println(hotspots);
    }
}
