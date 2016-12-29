package com.hartwig.healthchecker.common.checks;

import java.util.Optional;

import org.junit.Assert;
import org.junit.Test;

public class CheckTypeTest {

    @Test
    public void getByTypeSuccess() {
        final Optional<CheckType> checkType = CheckType.getByCategory("somatic");
        assert checkType.isPresent();
        Assert.assertTrue(checkType.get() == CheckType.SOMATIC);
    }

    @Test
    public void getByTypeFailures() {
        final Optional<CheckType> checkType = CheckType.getByCategory("does not exist");
        Assert.assertFalse(checkType.isPresent());
    }
}