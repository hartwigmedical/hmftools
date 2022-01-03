
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - Health Checker"

dependencies {
    implementation(project(":hmf-common"))
    testImplementation(libs.junit)
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
}

application {
    mainClass.set("com.hartwig.hmftools.healthchecker.HealthChecksApplication")
}
