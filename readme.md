# Métodos Numéricos Aplicados al Problema de los Tres Cuerpos

## Descripción

Este proyecto buscar resolver el clásico problema de los tres cuerpos en la mecánica celeste a través de técnicas de simulación numérica. Mediante el uso de varios métodos numericos.

## Instalación y Ejecución

### Configuración del Entorno

Se recomienda el uso de un entorno virtual para evitar conflictos con las bibliotecas existentes. Puedes crear y activar un entorno virtual siguiendo estos pasos:

1. **Creación del Entorno Virtual:**

    En el directorio de tu proyecto, ejecuta:

    ```bash
    python -m venv venv
    ```

2. **Activación del Entorno Virtual:**

    - En Windows, usa:

      ```bash
      venv\Scripts\activate
      ```

    - En Unix o MacOS, usa:

      ```bash
      source venv/bin/activate
      ```

3. **Navegación al Directorio del Proyecto:**

    Accede al directorio que contiene el código fuente:

    ```bash
    cd code
    ```

4. **Instalación de Dependencias:**

    Con el entorno virtual activado, instala las dependencias necesarias:

    ```bash
    pip install -r requirements.txt
    ```

### Ejecución del Programa

Una vez configurado y activado tu entorno virtual, ejecuta el módulo principal:

```bash
python main.py
```

La ejecución de este programa iniciará una serie de simulaciones para diversas órbitas, empleando distintos métodos numéricos y variando los intervalos de tiempo de paso.

## ODESolvers: Paquete de Métodos Numéricos

Este proyecto incluye un módulo `ODESolvers` que puede ser utilizado para resolver ecuaciones diferenciales ordinarias (ODEs) utilizando varios métodos de Runge-Kutta. Aquí hay un ejemplo de cómo usarlo para resolver la ODE `y' = y - t^2 + 1` en el intervalo de tiempo de `0` a `2` con un valor inicial de `0.5` y un tamaño de paso de tiempo de `0.2` utilizando el método de Euler:

```python
from ODESolvers.methods import rk_methods
from ODESolvers.rungekutta import RKMethod

def fun(t, y):
    return y - t**2 + 1

y0 = 0.5
rk = RKMethod(rk_methods["euler"], fun, 0, 2, [y0], 0.2)
a, b = rk.run()
```

Aquí está la explicación detallada:

- **Importación de Métodos de Runge-Kutta:**

  Importamos el diccionario `rk_methods` del módulo `ODESolvers.methods`. Este contiene varios métodos de Runge-Kutta, incluyendo el método de Euler.

- **Uso de la Clase RKMethod:**

  Importamos la clase `RKMethod` del módulo `ODESolvers.rungekutta`. Esta clase permite crear un objeto capaz de resolver una Ecuación Diferencial Ordinaria (ODE) utilizando un método de Runge-Kutta específico.

- **Definición de la Función de la ODE:**

  Definimos la función `fun(t, y)`, que representa la ODE a resolver. Por ejemplo, la ODE puede ser `y' = y - t^2 + 1`.

- **Establecimiento del Valor Inicial:**

  Fijamos el valor inicial `y0` de la ODE en 0.5.

- **Creación del Objeto RKMethod:**

  Creamos un objeto `rk` de la clase `RKMethod`. Le pasamos el método de Euler (`rk_methods["euler"]`), la función `fun`, el intervalo de tiempo de la simulación (de 0 a 2), el valor inicial `[y0]`, y el tamaño del paso de tiempo 0.2.

- **Ejecución del Método run():**

  Ejecutamos el método `run()` del objeto `rk` para resolver la ODE. Este método devuelve dos arreglos: el primero contiene los tiempos en los que se calculó la solución, y el segundo los valores de la solución en esos tiempos. Estos se almacenan en `a` y `b`, respectivamente.
