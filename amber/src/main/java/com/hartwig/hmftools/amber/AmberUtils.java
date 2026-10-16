package com.hartwig.hmftools.amber;

import static java.lang.Double.isFinite;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.AmberConstants.QUAL_FILTERED_THRESHOLD;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.Doubles;

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

    public static boolean aboveQualFilter(final PositionEvidence posEvidence)
    {
        if(posEvidence.ReadDepth == 0)
            return false;

        int filteredCount = posEvidence.BaseQualFiltered + posEvidence.MapQualFiltered + posEvidence.SeqTechFiltered;
        double filteredPercent = filteredCount / (double)posEvidence.ReadDepth;
        return filteredPercent < QUAL_FILTERED_THRESHOLD;
    }

    public static boolean aboveTumorOnlyMinSupport(final PositionEvidence posEvidence, final int minSupport)
    {
        return posEvidence.RefSupport >= minSupport && posEvidence.AltSupport >= minSupport;
    }

    public static boolean aboveNoiseFloor(final TumorBAF tumorBAF, final double noiseFloor)
    {
        if(tumorBAF.TumorEvidence.lohRecovered())
            return true;

        return (isFinite(tumorBAF.refFrequency()) && Doubles.greaterOrEqual(tumorBAF.refFrequency(), noiseFloor))
            && (isFinite(tumorBAF.altFrequency()) && Doubles.greaterOrEqual(tumorBAF.altFrequency(), noiseFloor));
    }

    public static AmberBAF fromTumorBaf(final TumorBAF tumorBAF)
    {
        int tumorAltCount = tumorBAF.TumorEvidence.AltSupport;
        double tumorBaf = tumorAltCount / (double) (tumorAltCount + tumorBAF.TumorEvidence.RefSupport);
        int normalAltCount = tumorBAF.NormalAltSupport;
        double normalBaf = normalAltCount / (double) (normalAltCount + tumorBAF.NormalRefSupport);
        boolean lohRecovered = tumorBAF.TumorEvidence.lohRecovered();

        return new AmberBAF(
                tumorBAF.chromosome(), tumorBAF.position(), tumorBaf, tumorBAF.TumorEvidence.ReadDepth, normalBaf,
                tumorBAF.NormalReadDepth, lohRecovered);
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
