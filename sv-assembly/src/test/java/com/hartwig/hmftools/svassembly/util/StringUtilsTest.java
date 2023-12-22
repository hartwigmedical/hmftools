package com.hartwig.hmftools.svassembly.util;

import org.junit.Test;

import static org.assertj.core.api.Assertions.assertThat;

import com.hartwig.hmftools.svassembly.util.StringUtils;

public class StringUtilsTest {
    @Test
    public void handlesNormalCases() {
        assertThat(StringUtils.toSnakeCase("ratherSimpleTest")).isEqualTo("rather_simple_test");
    }

    @Test
    public void handlesShortWords() {
        assertThat(StringUtils.toSnakeCase("thisIsATest")).isEqualTo("this_is_a_test");
    }

    @Test
    public void handlesLeadingAndTrailingAcronyms() {
        assertThat(StringUtils.toSnakeCase("BAMFileBAM")).isEqualTo("bam_file_bam");
    }

    @Test
    public void allUppercase() {
        assertThat(StringUtils.toSnakeCase("TODAY")).isEqualTo("today");
    }
}