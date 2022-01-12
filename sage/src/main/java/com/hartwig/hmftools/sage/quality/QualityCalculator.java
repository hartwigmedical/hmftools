package com.hartwig.hmftools.sage.quality;

import java.util.Map;

import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.config.QualityConfig;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import htsjdk.samtools.SAMRecord;

public class QualityCalculator
{
    private final QualityConfig mConfig;
    private final Map<String,QualityRecalibrationMap> mQualityRecalibrationMap;
    private final IndexedBases mRefBases;

    public QualityCalculator(
            final QualityConfig config, final Map<String, QualityRecalibrationMap> qualityRecalibrationMap, final IndexedBases refBases)
    {
        mConfig = config;
        mQualityRecalibrationMap = qualityRecalibrationMap;
        mRefBases = refBases;
    }

    public double calculateQualityScore(
            final ReadContextCounter readContextCounter, int readBaseIndex, final SAMRecord record, int numberOfEvents)
    {
        double baseQuality = baseQuality(readContextCounter, readBaseIndex, record);
        int distanceFromReadEdge = readDistanceFromEdge(readContextCounter, readBaseIndex, record);

        int mapQuality = record.getMappingQuality();
        boolean properPairFlag = record.getReadPairedFlag() && record.getProperPairFlag();
        int modifiedMapQuality = mConfig.modifiedMapQuality(readContextCounter.variant(), mapQuality, numberOfEvents, properPairFlag);
        double modifiedBaseQuality = mConfig.modifiedBaseQuality(baseQuality, distanceFromReadEdge);

        return Math.max(0, Math.min(modifiedMapQuality, modifiedBaseQuality));
    }

    private double baseQuality(final ReadContextCounter readContextCounter, int readBaseIndex, SAMRecord record)
    {
        return !readContextCounter.variant().isIndel()
                ? baseQuality(readContextCounter, readBaseIndex, record, readContextCounter.variant().ref().length())
                : readContextCounter.readContext().avgCentreQuality(readBaseIndex, record);
    }

    private double baseQuality(final ReadContextCounter readContextCounter, int startReadIndex, final SAMRecord record, int length)
    {
        int maxIndex = Math.min(startReadIndex + length, record.getBaseQualities().length) - 1;
        int maxLength = maxIndex - startReadIndex + 1;

        QualityRecalibrationMap qrMap = mQualityRecalibrationMap.get(readContextCounter.Sample);

        double quality = Integer.MAX_VALUE;
        for(int i = 0; i < maxLength; i++)
        {
            int refPosition = readContextCounter.position() + i;
            int readIndex = startReadIndex + i;
            byte rawQuality = record.getBaseQualities()[readIndex];
            byte[] trinucleotideContext = mRefBases.trinucleotideContext(refPosition);

            if(qrMap != null)
            {
                double recalibratedQuality = qrMap.quality(
                        (byte) readContextCounter.ref().charAt(i), (byte) readContextCounter.alt().charAt(i),
                        trinucleotideContext, rawQuality);

                quality = Math.min(quality, recalibratedQuality);
            }
        }

        return quality;
    }

    private int readDistanceFromEdge(final ReadContextCounter readContextCounter, int readIndex, final SAMRecord record)
    {
        int index = readContextCounter.readContext().readBasesPositionIndex();
        int leftIndex = readContextCounter.readContext().readBasesLeftCentreIndex();
        int rightIndex = readContextCounter.readContext().readBasesRightCentreIndex();

        int leftOffset = index - leftIndex;
        int rightOffset = rightIndex - index;

        int adjustedLeftIndex = readIndex - leftOffset;
        int adjustedRightIndex = readIndex + rightOffset;

        return Math.max(0, Math.min(adjustedLeftIndex, record.getReadBases().length - 1 - adjustedRightIndex));
    }

}
