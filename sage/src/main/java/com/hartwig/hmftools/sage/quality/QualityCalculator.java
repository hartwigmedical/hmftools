package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.MAX_MAP_QUALITY;
import static com.hartwig.hmftools.sage.SageConstants.READ_EDGE_PENALTY_0;
import static com.hartwig.hmftools.sage.SageConstants.READ_EDGE_PENALTY_1;
import static com.hartwig.hmftools.sage.bqr.BqrConfig.useReadType;
import static com.hartwig.hmftools.sage.bqr.BqrRegionReader.extractReadType;
import static com.hartwig.hmftools.sage.evidence.ArtefactContext.NOT_APPLICABLE_BASE_QUAL;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.sage.bqr.BqrRecordMap;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import htsjdk.samtools.SAMRecord;

public class QualityCalculator
{
    private final QualityConfig mConfig;
    private final BqrRecordMap mQualityRecalibrationMap;
    private final MsiJitterCalcs mMsiJitterCalcs;
    private final RefSequence mRefBases;
    private final boolean mUseReadType;
    private final SequencingType mSequencingType;
    private final UltimaQualCalculator mUltimaQualCalculator;

    private static final int MAX_HIGHLY_POLYMORPHIC_GENES_QUALITY = 10;

    public QualityCalculator(
            final SageConfig config, final BqrRecordMap qualityRecalibrationMap, final RefSequence refBases,
            final RefGenomeInterface refGenome, final MsiJitterCalcs msiJitterCalcs)
    {
        mConfig = config.Quality;
        mUseReadType = useReadType(config);
        mSequencingType = config.Sequencing.Type;
        mQualityRecalibrationMap = qualityRecalibrationMap;
        mMsiJitterCalcs = msiJitterCalcs;

        mRefBases = refBases;

        mUltimaQualCalculator = mSequencingType == SequencingType.ULTIMA ? new UltimaQualCalculator(refGenome) : null;
    }

    public boolean ultimaEnabled() { return mUltimaQualCalculator != null; }
    public MsiJitterCalcs msiJitterCalcs() { return mMsiJitterCalcs; }

    public UltimaQualModel createUltimaQualModel(final SimpleVariant variant)
    {
        return mUltimaQualCalculator != null ? mUltimaQualCalculator.buildContext(variant) : null;
    }

    public static int modifiedMapQuality(
            final QualityConfig config, final BasePosition position, int mapQuality, double readEvents, boolean isImproperPair)
    {
        if(config.isHighlyPolymorphic(position))
        {
            return min(MAX_HIGHLY_POLYMORPHIC_GENES_QUALITY, mapQuality - config.FixedPenalty);
        }

        int improperPairPenalty = isImproperPair ? config.ImproperPairPenalty : 0;
        int eventPenalty = (int)round(max(0, readEvents - 1) * config.ReadEventsPenalty);

        int modifiedMapQuality = mapQuality - config.FixedPenalty - improperPairPenalty - eventPenalty;

        return config.MapQualityRatioFactor > 0 ? min(MAX_MAP_QUALITY, modifiedMapQuality) : modifiedMapQuality;
    }

    public static class QualityScores
    {
        public final double RawBaseQuality;
        public final double RecalibratedBaseQuality;
        public final int ModifiedMapQuality;
        public final double ModifiedBaseQuality;
        public final double ModifiedQuality;

        public QualityScores(
                double rawBaseQuality, double recalibratedBaseQuality, int modifiedMapQuality,
                double modifiedBaseQuality, double modifiedQuality)
        {
            RawBaseQuality = rawBaseQuality;
            RecalibratedBaseQuality = recalibratedBaseQuality;
            ModifiedMapQuality = modifiedMapQuality;
            ModifiedBaseQuality = modifiedBaseQuality;
            ModifiedQuality = modifiedQuality;
        }
    }

    public QualityScores calculateQualityScores(
            final ReadContextCounter readContextCounter, int readBaseIndex, final SAMRecord record, double numberOfEvents, double rawBaseQuality)
    {
        double baseQuality;

        if(readContextCounter.isIndel() || readContextCounter.artefactContext() != null
        || (readContextCounter.ultimaQualModel() != null && rawBaseQuality != ULTIMA_MAX_QUAL))
        {
            baseQuality = rawBaseQuality;
        }
        else
        {
            baseQuality = recalibratedBaseQuality(readContextCounter, readBaseIndex, record, readContextCounter.variant().ref().length());
        }

        int mapQuality = record.getMappingQuality();
        boolean isImproperPair = isImproperPair(record);

        int modifiedMapQuality = modifiedMapQuality(mConfig, readContextCounter.variant(), mapQuality, numberOfEvents, isImproperPair);

        double modifiedBaseQuality = baseQuality - mConfig.BaseQualityFixedPenalty;

        int readEdgePenalty = readEdgeDistancePenalty(readContextCounter, readBaseIndex, record);
        modifiedBaseQuality = modifiedBaseQuality - readEdgePenalty;

        double modifiedQuality = max(0, min(modifiedMapQuality, modifiedBaseQuality));

        return new QualityScores(
                rawBaseQuality, baseQuality, max(0, modifiedMapQuality), max(0.0, modifiedBaseQuality), modifiedQuality);
    }

    public static boolean isImproperPair(final SAMRecord record) { return record.getReadPairedFlag() && !record.getProperPairFlag(); }

    public static double rawBaseQuality(final ReadContextCounter readContextCounter, int readIndex, final SAMRecord record)
    {
        if(readContextCounter.ultimaQualModel() != null)
            return readContextCounter.ultimaQualModel().calculateQual(record, readIndex);

        if(readContextCounter.artefactContext() != null)
        {
            byte adjustBaseQual = readContextCounter.artefactContext().findApplicableBaseQual(record, readIndex);
            if(adjustBaseQual != NOT_APPLICABLE_BASE_QUAL)
                return adjustBaseQual;
        }

        if(readContextCounter.isIndel())
        {
            double avgCoreQuality = readContextCounter.readContextMatcher().averageCoreQuality(record, readIndex);
            double msiIndelErrorRate = readContextCounter.qualCache().msiIndelErrorRate();
            return msiIndelErrorRate > 0 ? min(avgCoreQuality, msiIndelErrorRate) : avgCoreQuality;
        }

        if(readContextCounter.isSnv())
            return record.getBaseQualities()[readIndex];

        int varLength = readContextCounter.variant().ref().length();

        // take the minimum base qual for MNVs
        double minBaseQual = record.getBaseQualities()[readIndex];

        for(int i = readIndex + 1; i < readIndex + varLength; ++i)
        {
            minBaseQual = min(minBaseQual, record.getBaseQualities()[i]);
        }

        return minBaseQual;
    }

    private double recalibratedBaseQuality(
            final ReadContextCounter readContextCounter, int startReadIndex, final SAMRecord record, int length)
    {
        BqrReadType readType = mUseReadType ? extractReadType(record, mSequencingType) : BqrReadType.NONE;

        if(readContextCounter.isSnv())
        {
            // simplified version of the MNV case below
            byte rawQuality = record.getBaseQualities()[startReadIndex];
            return readContextCounter.qualCache().getQual(rawQuality, readType, 0);
        }

        // MNV case
        int maxIndex = min(startReadIndex + length, record.getBaseQualities().length) - 1;
        int maxLength = maxIndex - startReadIndex + 1;

        double quality = Integer.MAX_VALUE;
        for(int i = 0; i < maxLength; i++)
        {
            int readIndex = startReadIndex + i;
            byte rawQuality = record.getBaseQualities()[readIndex];

            double recalibratedQual = readContextCounter.qualCache().getQual(rawQuality, readType, i);
            quality = min(quality, recalibratedQual);
        }

        return quality;
    }

    public byte[] getTrinucleotideContext(int refPosition)
    {
        return mRefBases.containsPosition(refPosition) ? mRefBases.trinucleotideContext(refPosition) : null;
    }

    public double lookupRecalibrateQuality(final byte[] trinucleotideContext, byte altBase, byte rawQuality, final BqrReadType readType)
    {
        if(rawQuality == 0)
            return 0; // never adjust a zero qual up

        if(mQualityRecalibrationMap == null)
            return rawQuality;

        return mQualityRecalibrationMap.getQualityAdjustment(trinucleotideContext[1], altBase, trinucleotideContext, rawQuality, readType);
    }

    private int readEdgeDistancePenalty(final ReadContextCounter readContextCounter, int readIndex, final SAMRecord record)
    {
        int minDistance;

        if(readContextCounter.isIndel())
        {
            VariantReadContext readContext = readContextCounter.readContext();
            int lowerVarIndex = readIndex - (readContext.VarIndex - readContext.AltIndexLower);
            int upperVarIndex = readIndex + (readContext.AltIndexUpper - readContext.VarIndex);

            minDistance = max(0, min(lowerVarIndex, record.getReadBases().length - 1 - upperVarIndex));
        }
        else
        {
            minDistance = max(0, min(readIndex, record.getReadBases().length - 1 - readIndex));
        }

        if(minDistance <= 0)
            return READ_EDGE_PENALTY_0;
        else if(minDistance == 1)
            return READ_EDGE_PENALTY_1;
        else
            return 0;
    }
}
