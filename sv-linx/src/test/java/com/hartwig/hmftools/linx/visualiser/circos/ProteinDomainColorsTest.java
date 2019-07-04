package com.hartwig.hmftools.linx.visualiser.circos;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import java.awt.Color;
import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class ProteinDomainColorsTest
{

    private static final Color COLOR0 = new Color(204,61,61);
    private static final Color COLOR1 = new Color(204,77,61);
    private static final Color COLOR2 = new Color(204,93,61);
    private static final Color COLOR3 = new Color(204,109,61);
    private static final Color COLOR4 = new Color(204,124,61);
    private static final Color COLOR5 = new Color(204,140,61);
    private static final Color COLOR6 = new Color(204,156,61);
    private static final Color COLOR7 = new Color(204,172,61);
    private static final Color COLOR8 = new Color(204,188,61);
    private static final Color COLOR9 = new Color(204,203,61);

    @Test
    public void testFloating()
    {

        final Set<String> domains = Sets.newHashSet("Cadherin", "PD1", "PD2", "PD3");
        final ProteinDomainColors domainColors = new ProteinDomainColors(domains);

        assertNotEquals(domainColors.color("Cadherin"), domainColors.color("PD3"));
        assertEquals(domainColors.color("PD3"), domainColors.color("PD3"));

        System.out.println(domainColors);

    }

}
