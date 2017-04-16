# Thinfilm

source : https://arxiv.org/abs/1603.02720

## Transfer Matrix

    ( ->1 )    1  (1                    -r21)  ( ->2 )
    (     ) = --- (                         )  (     )
    ( <-1 )   t12 (r12     t12 t21 - r12 r21)  ( <-2 )

## Dielectric interface Transfer Matrix

The transfer matrix between to medium with indices n1 and n2 and propagation
angle th1 and th2 (n1 sin(th1) = n2 sin(th2)) is given by

               1        (n1 cos(th1) + n2 cos(th2)    n1 cos(th1) - n2 cos(th2))
    M_S = ------------- (                                                      )
          2 n1 cos(th1) (n1 cos(th1) - n2 cos(th2)    n1 cos(th1) + n2 cos(th2))

for the S polarisation. And by

               1        (n2 cos(th1) + n1 cos(th2)    n2 cos(th1) - n1 cos(th2))
    M_P = ------------- (                                                      )
          2 n1 cos(th1) (n2 cos(th1) - n1 cos(th2)    n2 cos(th1) + n1 cos(th2))

for the P polarisation. For the interface, the transfer matrix reduces to

         1  (1   r12)
    M = --- (       )
        t12 (r12   1)

## Decomposition of the Transfer Matrix

                1       (n1 cos(th1)    1) (1                        1)
    M_S = ------------- (                ) (                          )
          2 n1 cos(th1) (n1 cos(th1)   -1) (n2 cos(th2)   -n2 cos(th2))

for the S polarisation. And

                 1       (n1    cos(th1)) (cos(th2)   cos(th2))
    M_P  = ------------- (              ) (                   )
           2 n1 cos(th1) (n1   -cos(th1)) (n2              -n2)

for P polarisation
