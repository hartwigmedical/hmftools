import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

// https://stackoverflow.com/questions/66450310/how-can-i-customize-a-kotlincompile-task-with-a-gradle-kotlin-buildsrc-plugin

plugins {
    kotlin("jvm")
    // include the java-conventions one as it includes some configurations that we want to apply
    id("com.hartwig.java-conventions")
}

// set jvm version to the one specified in gradle.properties
if (!hasProperty("java.version"))
    logger.error("java.version property missing")

val javaVersion = property("java.version").toString()

tasks.withType<KotlinCompile> {
    kotlinOptions {
        jvmTarget = javaVersion
    }
}
