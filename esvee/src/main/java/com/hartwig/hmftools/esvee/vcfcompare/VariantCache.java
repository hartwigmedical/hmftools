package com.hartwig.hmftools.esvee.vcfcompare;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.vcf.VCFHeader;

public class VariantCache
{
    private final GenotypeIds mGenotypeIds;
    private final StructuralVariantFactory mSvFactory;

    private final List<Variant> mSvData;
    private final Map<String,List<Breakend>> mChromosomeBreakends;

    public VariantCache(final String vcfFilename)
    {
        mSvData = Lists.newArrayList();
        mChromosomeBreakends = Maps.newHashMap();

        if(vcfFilename == null)
        {
            mGenotypeIds = null;
            mSvFactory = null;
            return;
        }

        VcfFileReader vcfFileReader = new VcfFileReader(vcfFilename);
        VCFHeader vcfHeader = vcfFileReader.vcfHeader();
        List<String> vcfSampleNames = vcfHeader.getGenotypeSamples();

        String tumorId = "";
        String referenceId = null;

        if(vcfSampleNames.size() > 1)
        {
            tumorId = vcfSampleNames.get(1);
            referenceId = vcfSampleNames.get(0);
        }
        else
        {
            tumorId = vcfSampleNames.get(0);
        }

        mGenotypeIds = new GenotypeIds(referenceId != null ? 1 : -1, 0, referenceId, tumorId);

        mSvFactory = new StructuralVariantFactory(new CompoundFilter(false));
        mSvFactory.setGenotypeOrdinals(mGenotypeIds.ReferenceOrdinal, mGenotypeIds.TumorOrdinal);

        vcfFileReader.iterator().forEach(x -> processVariant(x));

        buildBreakendMap();

        SV_LOGGER.debug("loaded {} SVs from vcf({})", mSvData.size(), vcfFilename);
    }

    public boolean valid() { return mGenotypeIds != null; }
    public GenotypeIds genotypeIds() { return mGenotypeIds; }

    public List<Breakend> getChromosomeBreakends(final String chromosome)
    {
        return mChromosomeBreakends.containsKey(chromosome) ? mChromosomeBreakends.get(chromosome) : Collections.emptyList();
    }

    public List<Variant> getSvList() { return mSvData; }
    public Map<String,List<Breakend>> getBreakendMap() { return mChromosomeBreakends; }

    public int sglCount() { return (int)mSvData.stream().filter(x -> x.isSgl()).count(); }
    public int svCount() { return (int)mSvData.stream().filter(x -> !x.isSgl()).count(); }
    public int incompleteSVs() { return mSvFactory.unmatched().size(); }

    private void processVariant(final VariantContext variant)
    {
        if(!HumanChromosome.contains(variant.getContig()))
            return;

        boolean isSgl = StructuralVariantFactory.isSingleBreakend(variant);

        if(isSgl)
        {
            StructuralVariant sv = mSvFactory.createSingleBreakend(variant);
            addSvData(new Variant(sv, mGenotypeIds));
            return;
        }

        String mateId = StructuralVariantFactory.mateId(variant);

        if(mateId == null)
            return;

        int currentSvCount = mSvFactory.results().size();
        mSvFactory.addVariantContext(variant);

        // check if both breakends have now been encountered
        if(currentSvCount == mSvFactory.results().size())
            return;

        StructuralVariant sv = popLastSv(); // get and clear from storage

        if(sv == null)
            return;

        addSvData(new Variant(sv, mGenotypeIds));
    }

    private StructuralVariant popLastSv()
    {
        if(mSvFactory.results().isEmpty())
            return null;

        StructuralVariant sv = mSvFactory.results().get(0);
        mSvFactory.results().remove(0);

        return sv;
    }

    private void addSvData(final Variant var)
    {
        mSvData.add(var);
    }

    private void buildBreakendMap()
    {
        for(Variant var : mSvData)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                Breakend breakend = var.breakends()[se];

                if(breakend == null)
                    continue;

                List<Breakend> breakends = mChromosomeBreakends.get(breakend.Chromosome);

                if(breakends == null)
                {
                    breakends = Lists.newArrayList();
                    mChromosomeBreakends.put(breakend.Chromosome, breakends);
                }

                breakends.add(breakend);
            }
        }

        for(List<Breakend> breakends : mChromosomeBreakends.values())
        {
            Collections.sort(breakends, new BreakendPositionComparator());
        }
    }

    public static class BreakendPositionComparator implements Comparator<Breakend>
    {
        public int compare(final Breakend first, final Breakend second)
        {
            if(first.Position == second.Position)
            {
                if(first.Orient == second.Orient)
                    return 0;
                else
                    return first.Orient.isForward() ? -1 : 1;
            }
            else
            {
                return first.Position < second.Position ? -1 : 1;
            }
        }
    }

    public void clear()
    {
        mChromosomeBreakends.clear();
        mSvData.clear();
    }
}
