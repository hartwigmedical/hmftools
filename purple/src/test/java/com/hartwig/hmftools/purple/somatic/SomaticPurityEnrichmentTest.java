package com.hartwig.hmftools.purple.somatic;


import static com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory.CLNSIG;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.KATAEGIS_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_SOMATIC_LIKELIHOOD;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.REPORTABLE_TRANSCRIPTS;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.REPORTABLE_TRANSCRIPTS_DELIM;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.SomaticLikelihood;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.util.Strings;
import org.junit.Test;
import static junit.framework.TestCase.assertEquals;

import htsjdk.variant.variantcontext.VariantContext;
import junit.framework.TestCase;

public class SomaticPurityEnrichmentTest extends TestCase
{
    @Test
    public void testCalculateBiallelic()
    {
        /* 
        input values:
            CN = 2.00
            MACN = 0.975
            VCN = 2.44
            alleleReadCount = 1.2238
        
        Expected biallelic probability = 0.033
         */
        
        // MACN = (1 - BAF) * CN <=> BAF = 1 - MACN / CN
        double correspondingBAF = 1 - (0.975 / 2.00);
                
        PurpleCopyNumber returnedCnAndMacn = getCnAndMacn(2.00, correspondingBAF);

        assertEquals( 0.975, returnedCnAndMacn.minorAlleleCopyNumber(), 0.01);
        assertEquals(2.00, returnedCnAndMacn.averageTumorCopyNumber(), 0.01);
        
    }
    

    private static PurpleCopyNumber getCnAndMacn(final double copyNumber, final double averageActualBAF)
    {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome("chr1")
                .start(1)
                .end(1000)
                .averageTumorCopyNumber(copyNumber)
                .segmentStartSupport(SegmentSupport.NONE)
                .segmentEndSupport(SegmentSupport.NONE)
                .method(CopyNumberMethod.UNKNOWN)
                .bafCount(0)
                .depthWindowCount(1)
                .gcContent(0)
                .minStart(1)
                .maxStart(10)
                .averageObservedBAF(0.5)
                .averageActualBAF(averageActualBAF)
                .build();
    }
    
    
    
    
    
    
    

}