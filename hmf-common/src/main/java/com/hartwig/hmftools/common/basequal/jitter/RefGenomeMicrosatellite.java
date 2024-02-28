package com.hartwig.hmftools.common.basequal.jitter;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import com.google.common.collect.ArrayListMultimap;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.util.StringUtil;

public class RefGenomeMicrosatellite
{
    public static final Logger sLogger = LogManager.getLogger(RefGenomeMicrosatellite.class);

    public final ChrBaseRegion genomeRegion;
    public final byte[] unit;
    public final int numRepeat;
    public double mappability = Double.NaN;

    public RefGenomeMicrosatellite(final ChrBaseRegion genomeRegion, final byte[] unit)
    {
        this.genomeRegion = genomeRegion;
        this.unit = unit;
        this.numRepeat = genomeRegion.baseLength() / unit.length;
    }

    public RefGenomeMicrosatellite(final ChrBaseRegion genomeRegion, byte unit)
    {
        this(genomeRegion, new byte[] { unit });
    }

    public RefGenomeMicrosatellite(final String chromosome, int start, int end, final byte[] unit)
    {
        this(new ChrBaseRegion(chromosome, start, end), unit);
    }

    public RefGenomeMicrosatellite(final String chromosome, int start, int end, final byte unit)
    {
        this(new ChrBaseRegion(chromosome, start, end), unit);
    }

    public String chromosome()
    {
        return genomeRegion.chromosome();
    }

    public int referenceStart()
    {
        return genomeRegion.start();
    }

    public int referenceEnd()
    {
        return genomeRegion.end();
    }

    public int baseLength()
    {
        return genomeRegion.baseLength();
    }

    public String unitString()
    {
        return StringUtil.bytesToString(unit);
    }

    @Override
    public String toString()
    {
        return genomeRegion.toString() + ' ' + numRepeat + " x " + unitString();
    }

    // filter the microsatellites such that each type of (unit, length) is approximately the target count
    static List<RefGenomeMicrosatellite> filterMicrosatellites(List<RefGenomeMicrosatellite> inputList, int targetCountPerType)
    {
        sLogger.info("filtering {} microsatellite sites, target count per type = {}", inputList.size(), targetCountPerType);

        // first put them all into a multimap
        ArrayListMultimap<Pair<String, Integer>, RefGenomeMicrosatellite> unitLengthMicrosatelliteMap = ArrayListMultimap.create();

        // should be able to use a groupby method in guava
        for(RefGenomeMicrosatellite microsatellite : inputList)
        {
            Pair<String, Integer> k = Pair.of(microsatellite.unitString(), microsatellite.numRepeat);
            unitLengthMicrosatelliteMap.put(k, microsatellite);
        }

        // use same seed for now
        Random random = new Random(0);
        List<RefGenomeMicrosatellite> filteredList = new ArrayList<>();

        for(Pair<String, Integer> msType : unitLengthMicrosatelliteMap.keySet())
        {
            List<RefGenomeMicrosatellite> l = unitLengthMicrosatelliteMap.get(msType);
            double frac = ((double)targetCountPerType) / l.size();
            if(frac < 1.0)
            {
                l.stream().filter(o -> random.nextDouble() <= frac).forEach(filteredList::add);
            }
            else
            {
                filteredList.addAll(l);
            }
        }

        sLogger.info("filtered {} microsatellite sites down to {}", inputList.size(), filteredList.size());

        return filteredList;
    }
}
