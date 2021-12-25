
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - Health Checker"

dependencies {
    implementation(project(":hmf-common"))
    testImplementation(libs.junit)
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.healthchecker.HealthChecksApplication")
}
