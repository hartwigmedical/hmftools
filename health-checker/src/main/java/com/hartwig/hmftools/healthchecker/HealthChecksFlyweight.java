package com.hartwig.hmftools.healthchecker;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import com.hartwig.hmftools.healthchecker.runners.CheckType;
import com.hartwig.hmftools.healthchecker.runners.ErrorHandlingChecker;
import com.hartwig.hmftools.healthchecker.runners.HealthChecker;
import com.hartwig.hmftools.healthchecker.exception.NotFoundException;
import com.hartwig.hmftools.healthchecker.resource.ResourceWrapper;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.reflections.Reflections;

final class HealthChecksFlyweight {

    private static final Logger LOGGER = LogManager.getLogger(HealthChecksFlyweight.class);

    private static final Map<CheckType, HealthChecker> FLYWEIGHT = new HashMap<>();
    private static final HealthChecksFlyweight INSTANCE = new HealthChecksFlyweight();
    private static final Reflections BASE = new Reflections("com.hartwig.hmftools.healthchecker");

    private static final Set<Class<? extends ErrorHandlingChecker>> BASE_SET = BASE.getSubTypesOf(
            ErrorHandlingChecker.class);

    static {
        BASE_SET.forEach(checker -> {
            try {
                final HealthChecker checkerInstance = checker.newInstance();
                final ResourceWrapper resourceWrapper = checker.getAnnotation(ResourceWrapper.class);
                final CheckType checkType = resourceWrapper.type();

                FLYWEIGHT.put(checkType, checkerInstance);
            } catch (InstantiationException | IllegalAccessException e) {
                LOGGER.error(String.format("Error occurred when instantiating checker. Error -> %s", e.getMessage()));
            }
        });
    }

    private HealthChecksFlyweight() {
    }

    @NotNull
    static HealthChecksFlyweight getInstance() {
        return INSTANCE;
    }

    @NotNull
    HealthChecker getChecker(@NotNull final String type) throws NotFoundException {
        final Optional<CheckType> checkType = CheckType.getByCategory(type);
        if (!checkType.isPresent()) {
            throw new NotFoundException(String.format("Invalid CheckType informed %s", type));
        }
        return FLYWEIGHT.get(checkType.get());
    }

    @NotNull
    Collection<HealthChecker> getAllCheckers() {
        return FLYWEIGHT.values();
    }
}
