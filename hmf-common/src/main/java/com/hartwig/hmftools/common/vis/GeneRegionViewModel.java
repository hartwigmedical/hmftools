package com.hartwig.hmftools.common.vis;

import com.hartwig.hmftools.common.region.BaseRegion;

import org.apache.commons.lang3.NotImplementedException;

public interface GeneRegionViewModel
{
    BaseRegion region();
    GeneRegionViewModel updateRegion(final BaseRegion newRegion);
    GeneRegionViewModel updateBasesBeforeRightInsert(int newBasesBeforeRightInsert);
    int basesBeforeRightInsert();

    record AminoAcidViewModel(BaseRegion region, int aminoAcidPos, char ref, char alt, int basesBeforeRightInsert, boolean isStart)
            implements GeneRegionViewModel
    {
        public boolean matchesRef()
        {
            return ref == alt;
        }

        @Override
        public GeneRegionViewModel updateRegion(final BaseRegion newRegion)
        {
            return new AminoAcidViewModel(newRegion, aminoAcidPos, ref, alt, basesBeforeRightInsert, isStart);
        }

        @Override
        public GeneRegionViewModel updateBasesBeforeRightInsert(int newBasesBeforeRightInsert)
        {
            return new AminoAcidViewModel(region, aminoAcidPos, ref, alt, newBasesBeforeRightInsert, isStart);
        }
    }

    record NonCodingExonicRegionViewModel(BaseRegion region, int basesBeforeRightInsert) implements GeneRegionViewModel
    {
        @Override
        public GeneRegionViewModel updateRegion(final BaseRegion newRegion)
        {
            return new NonCodingExonicRegionViewModel(newRegion, basesBeforeRightInsert);
        }

        @Override
        public GeneRegionViewModel updateBasesBeforeRightInsert(int newBasesBeforeRightInsert)
        {
            return new NonCodingExonicRegionViewModel(region, newBasesBeforeRightInsert);
        }
    }

    record IntronicRegionViewModel(BaseRegion region, int basesBeforeRightInsert) implements GeneRegionViewModel
    {
        @Override
        public GeneRegionViewModel updateRegion(final BaseRegion newRegion)
        {
            return new IntronicRegionViewModel(newRegion, basesBeforeRightInsert);
        }

        @Override
        public GeneRegionViewModel updateBasesBeforeRightInsert(int newBasesBeforeRightInsert)
        {
            return new IntronicRegionViewModel(region, newBasesBeforeRightInsert);
        }
    }

    record DelViewModel(BaseRegion region) implements GeneRegionViewModel
    {
        @Override
        public GeneRegionViewModel updateRegion(final BaseRegion newRegion)
        {
            return new DelViewModel(newRegion);
        }

        @Override
        public GeneRegionViewModel updateBasesBeforeRightInsert(int newBasesBeforeRightInsert)
        {
            throw new NotImplementedException("Not implemented");
        }

        @Override
        public int basesBeforeRightInsert()
        {
            throw new NotImplementedException("Not implemented");
        }
    }
}
