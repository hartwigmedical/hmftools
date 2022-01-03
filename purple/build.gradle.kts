
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - PURPLE"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.mockito.core)
    testImplementation(libs.junit)
}

application {
    mainClass.set("com.hartwig.hmftools.purple.PurpleApplication")
}
