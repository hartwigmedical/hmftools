package com.hartwig.hmftools.dnds.calcs;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.dnds.DndsCommon.DN_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.dnds.SomaticVariant;

public class DndsMutation
{
    private final SomaticVariant mVariant;
    private final String mDndsImpact;
    private final Impact mImpact;

    private static final Set<Impact> EXCLUDE_HOTSPOT = Sets.newHashSet(Impact.SYNONYMOUS, Impact.UNKNOWN);
    private static final Set<Impact> EXCLUDE_BIALLELIC = Sets.newHashSet(Impact.MISSENSE, Impact.SYNONYMOUS, Impact.UNKNOWN);

    public DndsMutation(final SomaticVariant variant, final String dndsImpact)
    {
        mVariant = variant;
        mDndsImpact = dndsImpact;
        mImpact = determineImpact(variant, dndsImpact);
    }

    public boolean isBiallelic() { return mVariant.Biallelic && !EXCLUDE_BIALLELIC.contains(mImpact); }
    public boolean isHotspot() { return mVariant.Hotspot && !EXCLUDE_HOTSPOT.contains(mImpact); }
    public boolean isKnownOncoDriver() { return mVariant.Hotspot || (mImpact == Impact.INFRAME && mVariant.RepeatCount < 8); }
    public boolean isKnownTsgDriver() { return mVariant.Hotspot || mVariant.Biallelic; }

    private static Impact determineImpact(final SomaticVariant variant, final String dndsImpact)
    {
        if(dndsImpact.equals("no-SNV") || dndsImpact.equals("NA"))
        {
            if(variant.CanonicalCodingEffect == CodingEffect.MISSENSE)
                return Impact.INFRAME;
            else if(variant.CanonicalCodingEffect == CodingEffect.NONSENSE_OR_FRAMESHIFT)
                return Impact.FRAMESHIFT;
            else if(variant.CanonicalCodingEffect == CodingEffect.SYNONYMOUS)
                return Impact.SYNONYMOUS;

            if(variant.WorstCodingEffect == CodingEffect.MISSENSE)
                return Impact.INFRAME;
            else if(variant.WorstCodingEffect == CodingEffect.NONSENSE_OR_FRAMESHIFT)
                return Impact.FRAMESHIFT;
            else if(variant.WorstCodingEffect == CodingEffect.SYNONYMOUS)
                return Impact.SYNONYMOUS;

            return Impact.UNKNOWN;
        }

        if(dndsImpact.equals("Missense"))
            return Impact.MISSENSE;
        else if(dndsImpact.equals("Nonsense") || dndsImpact.equals("Stop_loss"))
            return Impact.NONSENSE;
        else if(dndsImpact.equals("Synonymous"))
            return Impact.SYNONYMOUS;
        else if(dndsImpact.equals("Essential_Splice"))
            return Impact.SPLICE;
        else
            return Impact.UNKNOWN;
    }

    /*

    class DndsMutationComparator(private val knownFunction: (DndsMutation) -> Boolean) : Comparator<DndsMutation> {

    override fun compare(o1: DndsMutation, o2: DndsMutation): Int {
        fun Boolean.toInt() = if (this) 1 else 0

        if (o1.isHotspot.xor(o2.isHotspot)) {
            return o2.isHotspot.toInt() - o1.isHotspot.toInt()
        }

        val o1Known = knownFunction(o1)
        val o2Known = knownFunction(o2)
        if (o1Known.xor(o2Known)) {
            return o2Known.toInt() - o1Known.toInt()
        }

        return o1.impact.ordinal - o2.impact.ordinal
    }
     */

    public static List<DndsMutation> readVariants(final String filename)
    {
        List<DndsMutation> mutations = Lists.newArrayList();

        try
        {
            final List<String> lines = Files.readAllLines(new File(filename).toPath());

            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
            lines.remove(0);

            int chrIndex = fieldsIndexMap.get("Chromosome");
            int posIndex = fieldsIndexMap.get("Position");
            int refIndex = fieldsIndexMap.get("Ref");
            int altIndex = fieldsIndexMap.get("Alt");
            int geneIndex = fieldsIndexMap.get("Gene");
            int biIndex = fieldsIndexMap.get("Biallelic");
            int hotspotIndex = fieldsIndexMap.get("Hotspot");
            int wceIndex = fieldsIndexMap.get("WorstCodingEffect");
            int ceIndex = fieldsIndexMap.get("CanonicalCodingEffect");
            int rcIndex = fieldsIndexMap.get("RepeatCount");
            int dndsImpactIndex = fieldsIndexMap.get("DndsImpact");

            for(String line : lines)
            {
                final String[] values = line.split(TSV_DELIM, -1);

                SomaticVariant variant = new SomaticVariant(
                        values[chrIndex], Integer.parseInt(values[posIndex]), values[refIndex], values[altIndex], values[geneIndex],
                        Boolean.parseBoolean(values[biIndex]), Boolean.parseBoolean(values[hotspotIndex]), CodingEffect.valueOf(values[wceIndex]),
                        CodingEffect.valueOf(values[ceIndex]), Integer.parseInt(values[rcIndex]));

                mutations.add(new DndsMutation(variant, values[dndsImpactIndex]));
            }
        }
        catch (IOException e)
        {
            DN_LOGGER.error("failed to read DNDS mutations: {}", e.toString());
            return null;
        }

        return mutations;
    }

}
