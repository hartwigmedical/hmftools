package com.hartwig.hmftools.svanalysis.analysis;

import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.SV_GROUP_ENCLOSED;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.SV_GROUP_ENCLOSING;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.SV_GROUP_NEIGHBOURS;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.SV_GROUP_OVERLAP;
import static com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator.NO_LINE_ELEMENT;

import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import com.hartwig.hmftools.svanalysis.types.SvFootprint;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvClusteringMethods {

    private static final Logger LOGGER = LogManager.getLogger(SvClusteringMethods.class);

    private final SvUtilities mUtils;

    private Map<String, Integer> mChrArmSvCount;
    private Map<String, Double> mChrArmSvExpected;
    private Map<String, Double> mChrArmSvRate;
    private double mMedianChrArmRate;

    private static final double REF_BASE_LENGTH = 10000000D;
    private static int MIN_FOOTPRINT_COUNT = 4;

    public SvClusteringMethods(final SvUtilities clusteringUtils)
    {
        mUtils = clusteringUtils;
        mChrArmSvCount = Maps.newHashMap();
        mChrArmSvExpected = Maps.newHashMap();
        mChrArmSvRate = Maps.newHashMap();
        mMedianChrArmRate = 0;
    }

    public void findLinkedSVsByDistance(SvCluster cluster, List<SvClusterData> unassignedVariants) {
        // look for any other SVs which form part of this cluster based on proximity
        int currentIndex = 0;

        while (currentIndex < unassignedVariants.size()) {
            SvClusterData currentVar = unassignedVariants.get(currentIndex);

            // compare with all other SVs in this cluster
            boolean matched = false;
            for (SvClusterData otherVar : cluster.getSVs())
            {
                // test each possible linkage
                if (!mUtils.areVariantsLinkedByDistance(currentVar, otherVar))
                {
                    //LOGGER.debug("non-linked SVs: v1({}) and v2({})", currentVar.posId(), otherVar.posId());
                    continue;
                }

                cluster.addVariant(currentVar);
                // LOGGER.debug("cluster({}) add matched variant({}), totalCount({})", cluster.getId(), currentVar.id(), cluster.getSVs().size());

                matched = true;
                break;
            }

            if(matched)
            {
                unassignedVariants.remove(currentIndex);

                // as soon as a new SV is added to this cluster, need to start checking from the beginning again
                currentIndex = 0;
            }
            else
            {
                ++currentIndex;
            }
        }
    }

    public void findFootprints(final String sampleId, SvCluster cluster)
    {
        if(cluster.getCount() < MIN_FOOTPRINT_COUNT * 2)
            return;

        List<SvClusterData> crossArmSvs = Lists.newArrayList();
        List<SvFootprint> footprints = Lists.newArrayList();

        int nextFootprintId = 1;

        for(SvClusterData var : cluster.getSVs())
        {
            // put to one side cross-chr, cross-arm and long-spanning SVs
            if(var.type() == StructuralVariantType.BND || !var.isLocal()
            || var.position(false) - var.position(true) > mUtils.getBaseDistance() * 10)
            {
                crossArmSvs.add(var);
                continue;
            }

            // attempt to find a footprint within proximity of this SV, otherwise make new one
            boolean found = false;

            for(SvFootprint footprint : footprints)
            {
                if (!mUtils.areAnyWithinRange(var.position(true), var.position(false), footprint.posStart(), footprint.posEnd(), mUtils.getBaseDistance())) {
                    continue;
                }

                found = true;
                footprint.addVariant(var, false);
                break;
            }

            if(!found)
            {
                SvFootprint footprint = new SvFootprint(nextFootprintId++, var.chromosome(true), var.getStartArm());
                footprint.addVariant(var, false);
                footprints.add(footprint);
            }
        }

        // now merge any footprints which are close to each other
        for(int i = 0; i < footprints.size(); ++i) {

            SvFootprint fp1 = footprints.get(i);

            for (int j = i + 1; j < footprints.size(); ) {
                SvFootprint fp2 = footprints.get(j);

                if (!mUtils.areAnyWithinRange(fp1.posStart(), fp1.posEnd(), fp2.posStart(), fp2.posEnd(), mUtils.getBaseDistance())) {
                    ++j;
                    continue;
                }

                LOGGER.debug("merge footprints {} and {}", fp1.posId(), fp2.posId());

                for (SvClusterData var : fp2.getSVs()) {
                    fp1.addVariant(var, false);
                }

                footprints.remove(j); // and keep index the same
            }
        }

        // remove any small footprints
        for(int i = 0; i < footprints.size();) {

            SvFootprint footprint = footprints.get(i);

            if (footprint.getCount(false) < MIN_FOOTPRINT_COUNT) {
                footprints.remove(i);
                continue;
            }

            ++i;
        }

        // finally assign spanning SVs and use these to bridge footprints if possible
        for(SvClusterData var : crossArmSvs)
        {
            SvFootprint startFootprint = null;
            SvFootprint endFootprint = null;

            for(SvFootprint footprint : footprints)
            {
                for(int a = 0; a < 2; ++a)
                {
                    boolean useStart = (a == 0);

                    if (!mUtils.areAnyWithinRange(var.position(useStart), var.position(useStart), footprint.posStart(), footprint.posEnd(), mUtils.getBaseDistance()))
                    {
                        continue;
                    }

                    if(useStart)
                        startFootprint = footprint;
                    else
                        endFootprint = footprint;

                    footprint.addVariant(var, true);
                }

                if(startFootprint != null && endFootprint != null)
                {
                    if(startFootprint == endFootprint) {
                        LOGGER.error("footprint{} matched both ends of spanningSV({})", startFootprint.posId(), var.posId());
                        continue;
                    }

                    // skip if already linked
                    if(startFootprint.getLinkedFootprints().contains(endFootprint))
                        continue;

                    LOGGER.debug("footprints {} and {} linked by spanningSV({})", startFootprint.posId(), endFootprint.posId(), var.posId());

                    // create links between these footprints
                    startFootprint.addLinkedFootprint(endFootprint);
                    endFootprint.addLinkedFootprint(startFootprint);
                }
            }
        }

        // log the final collection of FPs
        // log the remaining ones
        for(final SvFootprint footprint : footprints) {

            LOGGER.info("sample({}) cluster({}) footprint({}) has SVs({}) spanSVs({}) linkedFPs({})",
                    sampleId, cluster.getId(), footprint.posId(), footprint.getCount(false),
                    footprint.getSpanningSVs().size(), footprint.getLinkedFootprints().size());

            for (SvClusterData var : footprint.getSVs()) {

                LOGGER.info("footprint({}) has sv({})", footprint.posId(), var.posId());
            }

            for (SvClusterData var : footprint.getSpanningSVs()) {

                LOGGER.info("footprint({}) has spanningSV({})", footprint.posId(), var.posId());
            }

            for (SvFootprint lnkFP : footprint.getLinkedFootprints()) {

                LOGGER.info("footprint({}) linked to other fp({})", footprint.posId(), lnkFP.posId());
            }
        }

    }

    public void setChromosomalArmStats(final List<SvClusterData> allVariants)
    {
        mChrArmSvCount.clear();
        mChrArmSvExpected.clear();
        mChrArmSvRate.clear();
        mMedianChrArmRate = 0;

        // form a map of unique arm to SV count
        for(final SvClusterData var : allVariants)
        {
            String chrArmStart = mUtils.getVariantChrArm(var,true);
            String chrArmEnd = mUtils.getVariantChrArm(var,false);

            // ensure an entry exists
            if (!mChrArmSvCount.containsKey(chrArmStart)) {
                mChrArmSvCount.put(chrArmStart, 0);
            }

            // exclude LINE elements from back-ground rates
            if(var.isStartLineElement().equals(NO_LINE_ELEMENT)) {
                mChrArmSvCount.replace(chrArmStart, mChrArmSvCount.get(chrArmStart) + 1);
            }

            if(!var.isNullBreakend())
            {
                if (!chrArmStart.equals(chrArmEnd) && !mChrArmSvCount.containsKey(chrArmEnd))
                {
                    mChrArmSvCount.put(chrArmEnd, 0);
                }

                if (var.isEndLineElement().equals(NO_LINE_ELEMENT))
                {
                    mChrArmSvCount.replace(chrArmEnd, mChrArmSvCount.get(chrArmEnd) + 1);
                }
            }
        }

        // now determine the background rate by taking the median value from amongst the arms
        // factoring in the arms which have no Q (14-16, 21-22) and excluding the X & Ys
        for(Map.Entry<String, Integer> entry : mChrArmSvCount.entrySet()) {

            final String chrArm = entry.getKey();
            final String chromosome = mUtils.getChrFromChrArm(chrArm);
            final String arm = mUtils.getArmFromChrArm(chrArm);

            long chrArmLength = mUtils.getChromosomalArmLength(chromosome, arm);
            int svCount = entry.getValue();
            double ratePerLength = svCount / (chrArmLength / REF_BASE_LENGTH); // the factor isn't important

            mChrArmSvRate.put(chrArm, ratePerLength);
            // LOGGER.debug("chrArm({}) ratePerMill({}) from count({}) length({})", chrArm, ratePerLength, svCount, chrArmLength);
        }

        mChrArmSvRate = sortByValue(mChrArmSvRate, false);

        mMedianChrArmRate = 0;
        int chrArmIndex = 0;
        for(Map.Entry<String, Double> entry : mChrArmSvRate.entrySet())
        {
            // LOGGER.debug("chrArm({}: {}) svRate({})", chrArmIndex, entry.getKey(), entry.getValue());

            if(chrArmIndex == 20)
                mMedianChrArmRate = entry.getValue();

           ++chrArmIndex;
        }

        LOGGER.debug("median SV rate({})", mMedianChrArmRate);

        // now create another map of expected SV count per arm using the median rate
        for(Map.Entry<String, Double> entry : mChrArmSvRate.entrySet()) {

            final String chrArm = entry.getKey();
            final String chromosome = mUtils.getChrFromChrArm(chrArm);
            final String arm = mUtils.getArmFromChrArm(chrArm);

            long chrArmLength = mUtils.getChromosomalArmLength(chromosome, arm);
            double expectedSvCount = (int)Math.round((chrArmLength / REF_BASE_LENGTH) * mMedianChrArmRate);
            LOGGER.debug("chrArm({}) expectedSvCount({}) vs actual({})", chrArm, expectedSvCount, mChrArmSvCount.get(chrArm));

            mChrArmSvExpected.put(chrArm, expectedSvCount);
        }
    }

    public String getChrArmData(final SvClusterData var)
    {
        String chrArmStart = mUtils.getVariantChrArm(var,true);

        boolean hasEnd = !var.isNullBreakend();
        String chrArmEnd = hasEnd ? mUtils.getVariantChrArm(var,false) : "";

        // report Start SV count : Expected SV Count : End SV Count : Expected SV Count
        return String.format("%d,%.2f,%d,%.2f",
                mChrArmSvCount.get(chrArmStart), mChrArmSvExpected.get(chrArmStart),
                hasEnd ? mChrArmSvCount.get(chrArmEnd) : 0, hasEnd ? mChrArmSvExpected.get(chrArmEnd) : 0.0);
    }

    private static Map<String, Double> sortByValue(Map<String, Double> unsortMap, final boolean order)
    {
        List<Map.Entry<String, Double>> list = new LinkedList<Map.Entry<String, Double>>(unsortMap.entrySet());

        Collections.sort(list, new Comparator<Map.Entry<String, Double>>()
        {
            public int compare(Map.Entry<String, Double> o1,
                    Map.Entry<String, Double> o2)
            {
                if (order)
                {
                    return o1.getValue().compareTo(o2.getValue());
                }
                else
                {
                    return o2.getValue().compareTo(o1.getValue());

                }
            }
        });

        // Maintaining insertion order with the help of LinkedList
        Map<String, Double> sortedMap = new LinkedHashMap<String, Double>();
        for (Map.Entry<String, Double> entry : list)
        {
            sortedMap.put(entry.getKey(), entry.getValue());
        }

        return sortedMap;
    }

}
