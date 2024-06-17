package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class AmberUtils
{
    public static boolean isValid(final AmberBAF baf)
    {
        return Double.isFinite(baf.tumorBAF()) & Double.isFinite(baf.normalBAF());
    }

    public static List<ChrBaseRegion> loadBedFromResource(String resourcePath)
    {
        List<ChrBaseRegion> genomeRegions = Lists.newArrayList();

        java.io.InputStream bedStream = AmberUtils.class.getClassLoader().getResourceAsStream(resourcePath);

        if(bedStream == null)
        {
            AMB_LOGGER.error("unable to find resource bed file: {}", resourcePath);
            throw new RuntimeException("unable to find resource bed file: " + resourcePath);
        }

        List<String> fileContents = new BufferedReader(new InputStreamReader(
                AmberUtils.class.getClassLoader().getResourceAsStream(resourcePath))).lines().collect(Collectors.toList());

        for(String line : fileContents)
        {
            String[] values = line.split(TSV_DELIM, 3);
            genomeRegions.add(new ChrBaseRegion(values[0], Integer.parseInt(values[1]) + 1, Integer.parseInt(values[2])));
        }

        return genomeRegions;
    }

    public static AmberBAF fromTumorBaf(final TumorBAF tumor)
    {
        int tumorAltCount = tumor.TumorEvidence.AltSupport;
        double tumorBaf = tumorAltCount / (double) (tumorAltCount + tumor.TumorEvidence.RefSupport);
        int normalAltCount = tumor.NormalAltSupport;
        double normalBaf = normalAltCount / (double) (normalAltCount + tumor.NormalRefSupport);

        return new AmberBAF(tumor.chromosome(), tumor.position(), tumorBaf, tumor.TumorEvidence.ReadDepth, normalBaf, tumor.NormalReadDepth);
    }

    public static AmberBAF fromBaseDepth(final PositionEvidence baseDepth)
    {
        int normalAltCount = baseDepth.AltSupport;
        double normalBaf = normalAltCount / (double) (normalAltCount + baseDepth.RefSupport);

        return new AmberBAF(baseDepth.Chromosome, baseDepth.Position, -1, -1, normalBaf, baseDepth.ReadDepth);
    }

    public static AmberSite depthAsSite(final PositionEvidence baseDepth)
    {
        return new AmberSite(baseDepth.chromosome(), baseDepth.position(), baseDepth.ref(), baseDepth.alt(), false);
    }
}
