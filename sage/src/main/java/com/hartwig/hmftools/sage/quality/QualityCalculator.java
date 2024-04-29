package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.MAX_MAP_QUALITY;
import static com.hartwig.hmftools.sage.bqr.BqrConfig.useReadType;
import static com.hartwig.hmftools.sage.bqr.BqrRegionReader.extractReadType;
import static com.hartwig.hmftools.sage.evidence.ArtefactContext.NOT_APPLICABLE_BASE_QUAL;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.sage.bqr.BqrRecordMap;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import htsjdk.samtools.SAMRecord;

public class QualityCalculator
{
    private final QualityConfig mConfig;
    private final BqrRecordMap mQualityRecalibrationMap;
    private final RefSequence mRefBases;
    private final boolean mUseReadType;
    private final SequencingType mSequencingType;
    private final UltimaQualCalculator mUltimaQualCalculator;

    private static final int MAX_HIGHLY_POLYMORPHIC_GENES_QUALITY = 10;
    private static final String NO_MODEL_NAME = "";

    public QualityCalculator(
            final SageConfig config, final BqrRecordMap qualityRecalibrationMap, final RefSequence refBases,
            final RefGenomeInterface refGenome)
    {
        mConfig = config.Quality;
        mUseReadType = useReadType(config);
        mSequencingType = config.Sequencing.Type;
        mQualityRecalibrationMap = qualityRecalibrationMap;
        mRefBases = refBases;

        mUltimaQualCalculator = mSequencingType == SequencingType.ULTIMA ? new UltimaQualCalculator(refGenome) : null;
    }

    public UltimaQualModel createUltimateQualModel(final SimpleVariant variant)
    {
        return mUltimaQualCalculator != null ? mUltimaQualCalculator.buildContext(variant) : null;
    }

    public static int modifiedMapQuality(
            final QualityConfig config, final GenomePosition position, int mapQuality, double readEvents, boolean isImproperPair)
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

        if(mConfig.DistanceFromReadEdgeFactor > 0)
        {
            int distanceFromReadEdge = readDistanceFromEdge(readContextCounter, readBaseIndex, record);
            int readEdgePenalty = max(mConfig.DistanceFromReadEdgeFactor * distanceFromReadEdge - mConfig.DistanceFromReadEdgeFixedPenalty, 0);
            modifiedBaseQuality = min(modifiedBaseQuality, readEdgePenalty);
        }

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
            return readContextCounter.readContextMatcher().averageCoreQuality(record, readIndex);

        if(readContextCounter.isSnv())
            return record.getBaseQualities()[readIndex];

        int varLength = readContextCounter.variant().ref().length();

        double baseQualTotal = 0;

        for(int i = readIndex; i < readIndex + varLength; ++i)
        {
            baseQualTotal += record.getBaseQualities()[i];
        }

        return baseQualTotal / varLength;
    }

    private double recalibratedBaseQuality(
            final ReadContextCounter readContextCounter, int startReadIndex, final SAMRecord record, int length)
    {
        BqrReadType readType = mUseReadType ? extractReadType(record, mSequencingType) : BqrReadType.NONE;

        if(readContextCounter.isSnv())
        {
            // simplified version of the MNV case below
            byte rawQuality = record.getBaseQualities()[startReadIndex];
            return readContextCounter.bqrQualCache().getQual(rawQuality, readType, 0, this);
        }

        // MNV case
        int maxIndex = min(startReadIndex + length, record.getBaseQualities().length) - 1;
        int maxLength = maxIndex - startReadIndex + 1;

        double quality = Integer.MAX_VALUE;
        for(int i = 0; i < maxLength; i++)
        {
            int readIndex = startReadIndex + i;
            byte rawQuality = record.getBaseQualities()[readIndex];

            double recalibratedQual = readContextCounter.bqrQualCache().getQual(rawQuality, readType, i, this);
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

    private int readDistanceFromEdge(final ReadContextCounter readContextCounter, int readIndex, final SAMRecord record)
    {
        // calculate the left and right core positions in the context of this read
        int leftOffset = readContextCounter.readContext().leftCoreLength();
        int rightOffset = readContextCounter.readContext().rightCoreLength();

        int adjustedLeftIndex = readIndex - leftOffset;
        int adjustedRightIndex = readIndex + rightOffset;

        // take the smaller of the left and right core index
        return max(0, min(adjustedLeftIndex, record.getReadBases().length - 1 - adjustedRightIndex));
    }
}
