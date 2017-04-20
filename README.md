# Thinfilm

source : https://arxiv.org/abs/1603.02720

## Transfer Matrix

     / ->1 \     1   / 1                    -r21 \   / ->2 \
    (       ) = --- (                             ) (       )
     \ <-1 /    t12  \ r12     t12 t21 - r12 r21 /   \ <-2 /

## Dielectric interface Transfer Matrix

The transfer matrix between to medium with indices n1 and n2 and propagation
angle th1 and th2 (n1 sin(th1) = n2 sin(th2)) is given by

               1         /n1 cos(th1) + n2 cos(th2)    n1 cos(th1) - n2 cos(th2)\
    I_s = ------------- (                                                        )
          2 n1 cos(th1)  \n1 cos(th1) - n2 cos(th2)    n1 cos(th1) + n2 cos(th2)/

for the S polarisation. And by

               1         /n2 cos(th1) + n1 cos(th2)    n2 cos(th1) - n1 cos(th2)\
    I_p = ------------- (                                                        )
          2 n1 cos(th1)  \n2 cos(th1) - n1 cos(th2)    n2 cos(th1) + n1 cos(th2)/

for the P polarisation. For the interface, the transfer matrix reduces to

         1   /1   r12\
    I = --- (         )
        t12  \r12   1/

## Decomposition of the Transfer Matrix

                1        /n1 cos(th1)    1\   /1                       1\
    I_s = ------------- (                  ) (                           )
          2 n1 cos(th1)  \n1 cos(th1)   -1/   \n2 cos(th2)  -n2 cos(th2)/
          -------------- A_1s --------------   ----------- B_2s -----------

for the S polarisation. And

                1         /n1      cos(th1)\   /cos(th2)         cos(th2)\
    I_p  = ------------- (                  ) (                           )
           2 n1 cos(th1)  \n1     -cos(th1)/   \n2                    -n2/
           -------------- A_1p --------------   ----------- B_2s -----------

for P polarisation

## Propagation Transfer Matrix

Propagation in medium with complex refractive index. Where `n sin(th)` is real.

         / 1/psi     0 \
    J = (               )      with psi = exp(i 2pi n d / lambda_0 cos(th))
         \ 0       psi /

where `d` is the layer thickness and `lambda_0` is the wavelength.

## Multilayer Transfer Matrix

For a two layer :

    M =   I_01   J_1   I_12    J_2   I_23
        A_0 (B_1 J_1 A_1) (B_2 J_2 A_2) B_3  = A_0 C_1 C_2 B_3
            ---- C_1 ---- ---- C_2 ----

           /c    -i s / cos(th) / n \
    C_s = (                          )
           \-i s n cos(th)        c /
    
           /c      -i s cos(th) / n \
    C_p = (                          )
           \-i s n / cos(th)      c /
    
    with c = cos(2pi n d / lambda_0 cos(th))
         s = sin(2pi n d / lambda_0 cos(th))
