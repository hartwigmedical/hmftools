package com.hartwig.hmftools.linx.visualiser.circos;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import java.awt.Color;
import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class ProteinDomainColorsTest
{
    private static final Color COLOR0 = new Color(204, 61, 61);
    private static final Color COLOR1 = new Color(204, 90, 61);
    private static final Color COLOR2 = new Color(204, 118, 61);
    private static final Color COLOR3 = new Color(204, 147, 61);
    private static final Color COLOR4 = new Color(204, 175, 61);
    private static final Color COLOR5 = new Color(204, 204, 61);
    private static final Color COLOR6 = new Color(175, 204, 61);
    private static final Color COLOR7 = new Color(147, 204, 61);
    private static final Color COLOR8 = new Color(118, 204, 61);
    private static final Color COLOR9 = new Color(90, 204, 61);
    private static final Color COLOR10 = new Color(61, 204, 183);
    private static final Color COLOR11 = new Color(61, 75, 204);
    private static final Color COLOR12 = new Color(197, 61, 204);

    @Test
    public void testFloating()
    {

        final Set<String> domains = Sets.newHashSet("Cadherin", "PD1", "PD2", "PD3");
        final ProteinDomainColors domainColors = new ProteinDomainColors(domains);

        assertNotEquals(domainColors.color("Cadherin"), domainColors.color("PD3"));
        assertEquals(domainColors.color("PD3"), domainColors.color("PD3"));

    }

}
